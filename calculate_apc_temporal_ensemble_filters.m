function JOBFILE = calculate_apc_temporal_ensemble_filters(JOBFILE, PASS_NUMBER)

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

% Set the spectral weighting method to APC
spectral_weighting_method = 'apc';

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

% Read the APC parameters
apc_field = JOBFILE.Processing(PASS_NUMBER).Correlation.SpectralWeighting.APC;
if isfield(apc_field, 'Method')
    apc_method = lower(...
    JOBFILE.Processing(PASS_NUMBER).Correlation.SpectralWeighting.APC.Method);
else
    apc_method = 'magnitude';
end

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

% % Allocate arrays to hold vectors
% tx_temp = nan(num_regions_full, num_pairs_correlate);
% ty_temp = nan(num_regions_full, num_pairs_correlate);

% Allocate arrays to hold particle diameters
dp_y_full = nan(num_regions_full, num_pairs_correlate);
dp_x_full = nan(num_regions_full, num_pairs_correlate);

% Loop over all the images
for n = 1 : num_pairs_correlate
    
    % Image paths
    image_path_01 = JOBFILE.Processing(PASS_NUMBER).Frames.Paths{1}{n};
    image_path_02 = JOBFILE.Processing(PASS_NUMBER).Frames.Paths{2}{n};
    
    % Get image names
    [~, file_name_01] = fileparts(image_path_01);
    [~, file_name_02] = fileparts(image_path_02);
    
    % Inform the use
    fprintf(1, '%s Pass %d of %d, pair %d of %d\n', ...
        upper(spectral_weighting_method),  ...
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
        
        % Add the current complex correlation
        % to the ensemble complex correlation
        cross_corr_array(:, :, k) = ...
                    cross_corr_array(:, :, k) + ...
                    cross_corr_spectral;    
    end
    
    % End timer
    t2 = toc(t1);
    
    % Inform the user that correlations finished.
    fprintf(1, 'Correlated %d regions in %02f seconds.\n', ...
        num_regions_correlate, t2);
    
    % Save the correlation planes to the jobfile.
    % This is done to avoid passing multiple variables
    % to the different functions.
    % These will be deleted before the jobfile is saved.
    JOBFILE.Processing(PASS_NUMBER).Correlation.Planes = ...
        cross_corr_array;
end

% Measure the displacements from
% the correlation planes.
[~, ~, ...
    particle_diameters_y, ...
    particle_diameters_x] =...
    planes2vect(JOBFILE, PASS_NUMBER);

% Copy the particle diameters. 
% This is only done in case they need
% to be saved later if APC was done.
dp_y_full(grid_indices, :) = ...
    repmat(particle_diameters_x, [1, num_pairs_correlate]);
dp_x_full(grid_indices, :) = ...
    repmat(particle_diameters_y, [1, num_pairs_correlate]); 

% Delete the correlation planes
JOBFILE.Processing(PASS_NUMBER).Correlation = ...
    rmfield(JOBFILE.Processing(PASS_NUMBER).Correlation, 'Planes');
    
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