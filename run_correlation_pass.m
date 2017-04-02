function JOBFILE = run_correlation_pass(JOBFILE, PASS_NUMBER)

% Default to pass number one
if nargin < 2
    PASS_NUMBER = 1; 
end

% Create image path list.
JOBFILE = create_image_pair_path_list(JOBFILE, PASS_NUMBER);

% Grid the image
JOBFILE = grid_image(JOBFILE, PASS_NUMBER);

% Get the source field for any iterative methods
% This should run even if no iterative method was specified
JOBFILE = get_iterative_source_field(JOBFILE, PASS_NUMBER);

% Deform the grid if requested
% This will run even if DWO isn't specified
% % Note that right now, this won't do ANYTHING
% which is OK because I don't need DWO at this second.
JOBFILE = discrete_window_offset(JOBFILE, PASS_NUMBER);

% Allocate the results
JOBFILE = allocate_results(JOBFILE, PASS_NUMBER);

% Read the correlation method
% Read the spectral weighting method ('SCC', 'GCC', 'RPC', 'APC')
spetral_weighting_method = get_spectral_weighting_method(JOBFILE, PASS_NUMBER);

% Figure out if we're doing APC
isAPC = ~isempty(regexpi(spetral_weighting_method, 'apc'));

% Read the enemble parameters
%
% Ensemble domain string
ensemble_domain_string = lower(get_ensemble_domain(JOBFILE, PASS_NUMBER));
%
% Ensemble direction
% This specifies whether to do a temporal ensemble, 
% or a spatial ensemble, or no ensemble.
ensemble_direction_string = lower(get_ensemble_direction(JOBFILE, PASS_NUMBER));
%
% Parse the ensemble direction string to figure out
% which ensemble was specified
%
% Flag for "Don't do any ensemble"
do_no_ensemble = ~isempty(regexpi(lower(ensemble_direction_string), 'no'));
%
% Flag for "Do the temporal ensemble"
do_temporal_ensemble = or( ...
    ~isempty(regexpi(lower(ensemble_direction_string), 'tim')), ...
    ~isempty(regexpi(lower(ensemble_direction_string), 'tem')));
%
% Flag for "Do the spatial ensemble"
do_spatial_ensemble = ~isempty(regexpi(lower(ensemble_direction_string), 'spa'));

% Correlation grid points
grid_correlate_x = JOBFILE.Processing(PASS_NUMBER).Grid.Points.Correlate.X;
grid_correlate_y = JOBFILE.Processing(PASS_NUMBER).Grid.Points.Correlate.Y;

% Full list of grid points
grid_full_x = JOBFILE.Processing(PASS_NUMBER).Grid.Points.Full.X;
grid_full_y = JOBFILE.Processing(PASS_NUMBER).Grid.Points.Full.Y;

% Indices of grid points to correlate
grid_indices = JOBFILE.Processing(PASS_NUMBER). ...
    Grid.Points.Correlate.Indices;

% Number of pairs to correlate
num_regions_correlate = length(grid_correlate_x(:));

% Total number of grid points
num_regions_full = length(grid_full_x(:));

% Number of pairs to correlate
num_pairs_correlate = read_num_pairs(JOBFILE, PASS_NUMBER);

% Region sizes
[region_height, region_width] = get_region_size(JOBFILE, PASS_NUMBER);

% Determine if deform is requested
iterative_method = JOBFILE.Processing(1).Iterative.Method;

% Determine whether deform is requested
deform_requested = ~isempty(regexpi(iterative_method, 'def'));

% Subpixel fit parameters
%
% Estimated particle diameter
% % % TO DO: 
% % % Make a function called get_rpc_diameter(JOBFILE, PASS_NUMBER)
particle_diameter = ...
   JOBFILE.Processing(PASS_NUMBER). ...
   Correlation.EstimatedParticleDiameter;

% Make some filters. 
% These are to let the parfor loops run
rpc_filter = ...
    spectral_energy_filter(region_height, region_width, particle_diameter);

% GCC filter is just ones
gcc_filter = ones(region_height, region_width);

% Extract any method-specific parameters
% 
% Might want to move this
% out of the correlation pass code.
switch lower(spetral_weighting_method)
    case 'apc'
        apc_field = JOBFILE.Processing(PASS_NUMBER).Correlation.SpectralWeighting.APC;
        if isfield(apc_field, 'Method')
            apc_method = lower(...
            JOBFILE.Processing(PASS_NUMBER).Correlation.SpectralWeighting.APC.Method);
        else
            apc_method = 'magnitude';
        end       
end

% Fit method
% % % For now ignore this and always do three-point fit.
%%%%%%

% Make the spatial windows
[spatial_window_01, spatial_window_02] = ...
    make_spatial_windows(JOBFILE, PASS_NUMBER);

% Allocate the correlation planes
cross_corr_array = zeros(...
            region_height, region_width, num_regions_correlate);
        
% Count the number of passes
% that the job will perform
% (this is just for printing)
num_passes = determine_number_of_passes(JOBFILE);

% Allocate arrays to hold vectors
tx_temp = nan(num_regions_full, num_pairs_correlate);
ty_temp = nan(num_regions_full, num_pairs_correlate);

% Allocate arrays to hold particle diameters
dp_y_full = nan(num_regions_full, num_pairs_correlate);
dp_x_full = nan(num_regions_full, num_pairs_correlate);

% Read whether or not to do ensemble
do_ensemble = JOBFILE.Processing(PASS_NUMBER).Correlation.Ensemble.DoEnsemble;

% Calculate the temporal ensemble planes 

% Loop over all the images
for n = 1 : num_pairs_correlate
    
    % Unless the temporal ensemble was specified,
    % re-zero the arrays for holding the cross 
    % correlation planes. 
    switch lower(ensemble_direction_string)
        case {'spatial', 'none'}
            cross_corr_array(:) = 0;
        otherwise
            if not(do_ensemble)
                cross_corr_array(:) = 0;
            end
    end
    
    % Image paths
    image_path_01 = JOBFILE.Processing(PASS_NUMBER).Frames.Paths{1}{n};
    image_path_02 = JOBFILE.Processing(PASS_NUMBER).Frames.Paths{2}{n};
    
    % Get image names
    [~, file_name_01] = fileparts(image_path_01);
    [~, file_name_02] = fileparts(image_path_02);
    
    % Inform the use
    fprintf(1, '%s Pass %d of %d, pair %d of %d\n', ...
        upper(spetral_weighting_method),  ...
        PASS_NUMBER, num_passes, n, num_pairs_correlate);
    fprintf(1, '%s and %s\n', file_name_01, file_name_02);
   
    % Load the images
    image_raw_01 = double(imread(image_path_01));
    image_raw_02 = double(imread(image_path_02));

    % Before any deform action happens,
    % set the images to be processed 
    % (image_01 and image_02) to be equal
    % to the raw input images.
    % This is the "default" configuration.
    % If deform proceeds
    % then the data stored in
    % image_01 and image_02 (the raw images) will 
    % be replaced with the deformed images.
    image_01 = image_raw_01;
    image_02 = image_raw_02;

    % Initialize the source grid
    source_grid_x = JOBFILE.Processing(PASS_NUMBER).Grid.Points.Full.X;
    source_grid_y = JOBFILE.Processing(PASS_NUMBER).Grid.Points.Full.Y;
    source_displacement_x = nan(size(source_grid_x));
    source_displacement_y = nan(size(source_grid_y));
    
    % If doing deformation, 
    % then execute deformation
    if deform_requested
    
        % Read the deform parameters
        %
        % This is the grid from the previous pass, 
        % which will inform the deform method.
        source_grid_x = JOBFILE.Processing(PASS_NUMBER).Iterative. ...
            Source.Grid.X;
        source_grid_y = JOBFILE.Processing(PASS_NUMBER). ...
            Iterative.Source.Grid.Y;
        %
        % These are the displacements from the
        % previous pass, which will inform
        % the deform method.
        source_displacement_x = JOBFILE.Processing(PASS_NUMBER). ...
            Iterative.Source.Displacement.X(:, n);
        source_displacement_y = JOBFILE.Processing(PASS_NUMBER). ...
            Iterative.Source.Displacement.Y(:, n);

        % Deform method
        deform_interpolation_method = JOBFILE.Processing(1). ...
            Iterative.Deform.Interpolation;

        % Determine if the data present
        % would result in any deformation happening
        %
        % If all of the source displacements are zero,
        % then the images shouldn't be deformed, 
        % even if deformation is requested.
        % I (Matt Giarra) wrote the deform code
        % to check internally check for this condition
        % and to skip deformation if all the source
        % displacements are zero. Because of that,
        % this if-else block (the one you're reading right now)
        % is kind of redundant.
        % The reason I've put it here anyway
        % is to make the code more understandable:
        % I don't want it to look like deform is
        % always run no matter what. Without this block
        % that's how the code would read (I think), and you'd have
        % to go into the deform code itself to realize that 
        % zero-everywhere displacement fields result
        % in the image deformation getting skipped.
        deform_data_exist = or(any(source_displacement_x ~= 0), ...
            any(source_displacement_y ~= 0));

        % Determine whether to do deform
        deform_can_proceed = and(deform_requested, deform_data_exist);

        % Do the deform if conditions are satisfied
        if deform_can_proceed
            % Inform the user that
            % deform is happening.
            fprintf(1, 'Deforming images...\n');

            % Deform the images if requested.
            [image_01, image_02] = ...
                deform_image_pair(image_raw_01, image_raw_02, ...
                source_grid_x, source_grid_y, ...
                source_displacement_x, source_displacement_y, ...
                deform_interpolation_method);
        end
    end
    
    % Extract the interrogation regions
    % from the images. These are the regions
    % that will be cross-correlated
    % to provide the PIV displacement or 
    % velocity estimate.
    region_mat_01 = extract_sub_regions(image_01, ...
        [region_height, region_width], ...
        grid_correlate_x, grid_correlate_y);
    region_mat_02 = extract_sub_regions(image_02, ...
        [region_height, region_width], ...
        grid_correlate_x, grid_correlate_y);
    
    % Inform the user that correlations 
    % are about to happen.
    fprintf(1, 'Correlating image pair...\n');
    
    % Tic tock
    t1 = tic;
    
    % Loop over the regions
    for k = 1 : num_regions_correlate
        
        % Extract the subregions
        region_01 = region_mat_01(:, :, k);
        region_02 = region_mat_02(:, :, k);
       
        % Correlate the interrogation regions
        %
        % Take the Fourier transform of the first interroagion region
        FT_01 = fft2(spatial_window_01 .* ...
            (region_01 - mean(region_01(:))));
        %
        % Take the Fourier Transform of the second interroagion region
        FT_02 = fft2(spatial_window_02 .* ...
            (region_02 - mean(region_02(:))));
        
        % Calculate the cross correlation by
        % conjugate-multiplying the Fourier 
        % transforms of the two interrogation regions.
        cross_corr_spectral = fftshift((FT_01 .* conj(FT_02)));
        
        % Once the cross correlation has been calculated,
        % we must decide where to put it. This decision
        % changes depending on how and whether
        % the ensemble correlation was specified.
        % 
        % Here is the "decision tree" for where
        % to put the cross correlation plane. 
        % I'm making this up as I go so please be gentle.
        %
        % No ensemble:
        %   The array of correlation planes
        %   can be [m x n x num_grid_points]
        %   in size, and the vectors are
        %   calculated after that array has
        %   been populated.
        %
        % Temporal ensemble:
        %   The array of correlation planes
        %   can be [ m x n x num_grid_points]
        %   in size, and the vectors are
        %   calculated after all of the images
        %   have contributed to the correlation.
        %
        % Spatial ensemble:
        %   The array of correlation planes
        %   can be [m x n x num_grid_points]
        %   in size, and are all added together
        %   after the array has been populated
        %   and a single vector calculated 
        %   from the result. 
        %   A more efficient way to do this, rather
        %   than saving all of the planes to an array,
        %   would be to add all the planes together
        %   as they are calculated. However, I'm 
        %   choosing to save the planes to a larger
        %   array because doing so is compatible with 
        %   the other two ensemble options 
        %   ("none" and "temporal"). Most systems
        %   The memory requirements won't change from
        %   the other methods, and we don't have
        %   to make any decisions about how to 
        %   store the correlation planes internally.
          
        % Switch between spatial and spectral ensemble
        switch lower(ensemble_domain_string)
            
            case 'spectral'                   
                % For spectral ensemble, add the current complex
                % correlation to the ensemble complex correlation
                cross_corr_array(:, :, k) = ...
                    cross_corr_array(:, :, k) + ...
                    cross_corr_spectral;
                       
            otherwise
                % For spatial ensemble or no ensemble, take the inverse
                % FT of the spectral correlation (i.e., the spatial
                % correlation) and add this to the purely real 
                % ensemble correlation
                
                % Split the correlation into 
                % magnitude and phase. This will
                % be used for filtering 
                % (note this happens even if "no filtering"
                % is specified). 
                % Spectral correlation phase and magnitude
                [spectral_corr_phase, spectral_corr_mag] = ...
                    split_complex(cross_corr_spectral);
 
                % Switch between correlation methods
                switch lower(spetral_weighting_method)
                    case 'scc'                          
                        % For SCC, take the original magnitude
                        % as the spectral filter. 
                        spectral_filter = spectral_corr_mag;
                        
                    case 'apc'                    
                        % Automatically calculate the APC filter.
                        spectral_filter = ...
                            calculate_apc_filter(cross_corr_spectral, ...
                            particle_diameter, apc_method);
                        
                    case 'rpc'
                        % If RPC was specified then 
                        % set the spectral filter to 
                        % be the RPC filter.
                        spectral_filter = rpc_filter;
                        
                    case 'gcc'
                        % If GCC was specified then
                        % set the spectral filter to 
                        % be the GCC filter, which is
                        % ones everywhere (i.e., 
                        % no filter). Doing this here
                        % lets us always "filter" 
                        % the correlation later in the 
                        % code by multiplying by a "filter"
                        % array, without having to switch
                        % between "filtering methods"
                        spectral_filter = gcc_filter;
                end
                
                % Apply the phase filter to the spectral plane.
                % Note that this operation is legitimate even
                % if SCC or GCC is selected. The reason for this
                % is that if SCC is selected, then the "filter"
                % is just the original magnitude. 
                % I (Matt Giarra) have coded it this way to
                % simplify the control flow. Not sure if this will
                % turn out to be a nice way to do it.
                cross_corr_spectral_filtered = ...
                    spectral_filter .* spectral_corr_phase;
                
                % Take the inverse FT of the "filtered" correlation.
                cross_corr_spatial = fftshift(abs(ifft2(fftshift(...
                            cross_corr_spectral_filtered))));
                        
                % Add this correlation to the spatial ensemble
                cross_corr_array(:, :, k) = ...
                    cross_corr_array(:, :, k) + cross_corr_spatial;
        end 
    end
    
    % End timer on the loop over
    % number of regions to correlate.
    t2 = toc(t1);
    
    % Inform the user
    fprintf(1, 'Correlated %d regions in %02f seconds.\n', ...
        num_regions_correlate, t2);
    
    % Save the correlation planes to the jobfile.
    % This is done to avoid passing multiple variables
    % to the different functions.
    % These will be deleted before the jobfile is saved.
    JOBFILE.Processing(PASS_NUMBER).Correlation.Planes = ...
        cross_corr_array;
    
    % Extract displacements 
    % if the temporal enemble wasn't specified.
    if not(do_temporal_ensemble)
        % Measure the displacements from
        % the correlation planes.
        [ty, tx, ...
            particle_diameters_y, ...
            particle_diameters_x] = ...
            planes2vect(JOBFILE, PASS_NUMBER);
        
        % Add the measured displacements to the 
        % temporary array of displacements.
        ty_temp(grid_indices, n) = ty;
        tx_temp(grid_indices, n) = tx;
        
        % Add the particle diameters to the arrays
        % that hold particle diameters. This is only 
        % done in case the particle diameters
        % need to be saved for APC later. 
        % There's probably a better way to do this.
        dp_y_full(grid_indices, n) = particle_diameters_y;
        dp_x_full(grid_indices, n) = particle_diameters_x;       
    end
    
    % Print a carriage return after
    % the image pair is done processing
    fprintf(1, '\n');
end

% Extract the displacements if the 
% temporal correlation WAS specified.
if do_temporal_ensemble
    
    % Measure the displacements from
    % the correlation planes.
    [ty, tx, ...
        particle_diameters_y, ...
        particle_diameters_x] =...
        planes2vect(JOBFILE, PASS_NUMBER);
    
    % Add the measured displacements to the 
    % temporary array of displacements.
    % Since temporal ensemble collapses
    % all of the displacements 
    % onto a single grid, 
    % the vector fields for all
    % N of the image pairs can be said
    % to have the same displacement.
    % We do it this way so that 
    % subsequent passes can use
    % the results of this pass
    % even if the ensemble methods
    % are different.
    ty_temp(grid_indices, :) = repmat(ty, [1, num_pairs_correlate]);
    tx_temp(grid_indices, :) = repmat(tx, [1, num_pairs_correlate]);
    
    % Copy the particle diameters. 
    % This is only done in case they need
    % to be saved later if APC was done.
    dp_y_full(grid_indices, :) = ...
        repmat(particle_diameters_x, [1, num_pairs_correlate]);
    dp_x_full(grid_indices, :) = ...
        repmat(particle_diameters_y, [1, num_pairs_correlate]); 
end

% Allocate the "output" vectors
tx_full_output = nan(num_regions_full, num_pairs_correlate);
ty_full_output = nan(num_regions_full, num_pairs_correlate);

% Add source dispalacements 
% to the measured displacements.
for n = 1 : num_pairs_correlate
     
    % Read the source displacements
    source_displacement_x = JOBFILE.Processing(PASS_NUMBER). ...
        Iterative.Source.Displacement.X(:, n);
    source_displacement_y = JOBFILE.Processing(PASS_NUMBER). ...
        Iterative.Source.Displacement.Y(:, n);
    
    % Resample the source displacement
    % from the source grid
    % onto the current grid.
    [source_field_interp_tx, source_field_interp_ty] = ...
        resample_vector_field(...
            source_grid_x, source_grid_y, ...
            source_displacement_x, source_displacement_y, ...
            grid_full_x, grid_full_y);
        
    % Add the source displacement from
    % the iterative method to
    % the measured displacement
    tx_raw_full = tx_temp(:, n) + source_field_interp_tx;
    ty_raw_full = ty_temp(:, n) + source_field_interp_ty;
    
    % Save the current fields as 
    % the "output" fields.
    % This variable gets updated
    % after each post-processing step
    % is performed. The reason for this
    % is that the deform methods
    % should chose the final field
    % from the previous pass, and 
    % this way the "final" field
    % is always saved.
    tx_full_output(:, n) = tx_raw_full;
    ty_full_output(:, n) = ty_raw_full;

    % Save the results to the structure
    JOBFILE.Processing(PASS_NUMBER).Results.Displacement.Raw.X(:, n) = ...
        tx_raw_full;
    JOBFILE.Processing(PASS_NUMBER).Results.Displacement.Raw.Y(:, n) = ...
        ty_raw_full;   
end

% Do the validation if requested
do_validation = JOBFILE.Processing(PASS_NUMBER).Validation.DoValidation;

% Do validation if requested
if do_validation == true;
    
    % Inform the user
    fprintf(1, 'Validating vector fields...\n');
    
    % Determine whether all of the fields
    % need to be validated. If temporal ensemble
    % was done, then the different time steps
    % saved in the jobfile are all identical,
    % and only one validation needs to be performed.
    switch ensemble_direction_string
        case 'temporal'
            num_pairs_validate = 1;
        otherwise
            num_pairs_validate = num_pairs_correlate;    
    end
    
    % Loop over all the pairs.
    for n = 1 : num_pairs_validate
        
        % Read the raw displacements
        tx_raw_full = JOBFILE.Processing(PASS_NUMBER). ...
            Results.Displacement.Raw.X(:, n);
        ty_raw_full = JOBFILE.Processing(PASS_NUMBER). ...
            Results.Displacement.Raw.Y(:, n);
        
        % Calculate the validated field.
        [tx_val_full, ty_val_full, is_outlier_full] = ...
            validateField_prana(grid_full_x, grid_full_y, tx_raw_full, ty_raw_full);
        
         % Add validated vectors to the jobfile
        JOBFILE.Processing(PASS_NUMBER).Results.Displacement. ...
            Validated.X(:, n) = tx_val_full;
        JOBFILE.Processing(PASS_NUMBER).Results.Displacement. ...
            Validated.Y(:, n) = ty_val_full;
        JOBFILE.Processing(PASS_NUMBER).Results.Displacement. ...
            Validated.IsOutlier(:, n) = is_outlier_full;

        % Update the "output" fields
        tx_full_output(:, n) = tx_val_full;
        ty_full_output(:, n) = ty_val_full; 
    end
    
    % If fewer fields were validated
    % than correlated, that means
    % temporal ensemble was done.
    % In this case, copy the vectors
    % from the first frame to all
    % the other frames.
    if num_pairs_validate < num_pairs_correlate
        
        % Copy the horizontal displacements.
        JOBFILE.Processing(PASS_NUMBER).Results.Displacement. ...
            Validated.X = ...
            repmat(tx_val_full, [1, num_pairs_correlate]);
        
        % Copy the vertical displacements
        JOBFILE.Processing(PASS_NUMBER).Results.Displacement. ...
            Validated.Y = ...
            repmat(ty_val_full, [1, num_pairs_correlate]);
        
        % Copy the outlier flags
        JOBFILE.Processing(PASS_NUMBER).Results.Displacement. ...
            Validated.IsOutlier = ...
            repmat(is_outlier_full, [1, num_pairs_correlate]);
        
        % Copy the output vectors
        tx_full_output = repmat(tx_val_full, [1, num_pairs_correlate]);
        ty_full_output = repmat(ty_val_full, [1, num_pairs_correlate]);   
    end
end

% Determine whether to do smoothing.
do_smoothing = JOBFILE.Processing(PASS_NUMBER).Smoothing.DoSmoothing;

% Smoothing
if do_smoothing == true
    
    % Extract the smoothing parameters
    smoothing_kernel_diameter = JOBFILE.Processing(PASS_NUMBER).Smoothing.KernelDiameter;
    smoothing_kernel_std = JOBFILE.Processing(PASS_NUMBER).Smoothing.KernelStdDev;
    
    % Inform the user
    fprintf(1, 'Smoothing vector fields...\n');
    
    % Determine whether all of the fields
    % need to be smoothed. If temporal ensemble
    % was done, then the different time steps
    % saved in the jobfile are all identical,
    % and only one validation needs to be performed.
    switch ensemble_direction_string
        case 'temporal'
            num_pairs_smooth = 1;
        otherwise
            num_pairs_smooth = num_pairs_correlate;    
    end
    
    % Loop over all the pairs.
    for n = 1 : num_pairs_smooth
        tx_full_input = tx_full_output(:, n);
        ty_full_input = ty_full_output(:, n);
        
        % Calculate smoothed field
        tx_smoothed_full = smoothField(tx_full_input, smoothing_kernel_diameter, smoothing_kernel_std);
        ty_smoothed_full = smoothField(ty_full_input, smoothing_kernel_diameter, smoothing_kernel_std);
        JOBFILE.Processing(PASS_NUMBER).Results.Displacement.Smoothed.X(:, n) = tx_smoothed_full;
        JOBFILE.Processing(PASS_NUMBER).Results.Displacement.Smoothed.Y(:, n) = ty_smoothed_full;

        % Update the "output" fields
        tx_full_output(:, n) = tx_smoothed_full;
        ty_full_output(:, n) = ty_smoothed_full;
        
    end
    
    % If fewer fields were smoothed
    % than correlated, that means
    % temporal ensemble was done.
    % In this case, copy the vectors
    % from the first frame to all
    % the other frames.
    if num_pairs_smooth < num_pairs_correlate
        
        % Copy the horizontal displacements.
        JOBFILE.Processing(PASS_NUMBER).Results.Displacement. ...
            Smoothed.X = ...
            repmat(tx_smoothed_full, [1, num_pairs_correlate]);
        
        % Copy the vertical displacements
        JOBFILE.Processing(PASS_NUMBER).Results.Displacement. ...
            Smoothed.Y = ...
            repmat(ty_smoothed_full, [1, num_pairs_correlate]);
        
        % Copy the output vectors
        tx_full_output = repmat(tx_smoothed_full, [1, num_pairs_correlate]);
        ty_full_output = repmat(ty_smoothed_full, [1, num_pairs_correlate]);   
    end    
end

% Save the "output" fields to the jobfile results field.
JOBFILE.Processing(PASS_NUMBER).Results.Displacement.Final.X = tx_full_output;
JOBFILE.Processing(PASS_NUMBER).Results.Displacement.Final.Y = ty_full_output;

% Delete the correlation planes
JOBFILE.Processing(PASS_NUMBER).Correlation = ...
    rmfield(JOBFILE.Processing(PASS_NUMBER).Correlation, 'Planes');

% If APC was done then save the effective
% particle diameters to the pass results.
if isAPC
    
    % Add the horizontal (columns direction)
    % APC diameter to the results.
    JOBFILE.Processing(PASS_NUMBER).Results. ...
        Filtering.APC.Diameter.X = dp_x_full;
    
    % Add the vertical (rows direction)
    % APC diameter to the results.
    JOBFILE.Processing(PASS_NUMBER).Results. ...
        Filtering.APC.Diameter.Y = dp_y_full;
end


end













