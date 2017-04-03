function [TY, TX, PARTICLE_DIAMETER_Y, PARTICLE_DIAMETER_X] = planes2vect(JOBFILE, PASS_NUMBER)
% This function measures displacements (TY, TX) from
% a list of either real or complex correlation planes. 

% Read the ensemble direction
ensemble_direction_string = lower( ...
    get_ensemble_direction(JOBFILE, PASS_NUMBER));

% Extract the correlation planes from the job file
cross_corr_array = JOBFILE.Processing(PASS_NUMBER).Correlation.Planes;

% Check if all the imaginary components
% of the correlation array are equal to zero.
% If this is the case, then spectral plane fitting
% won't work.
planes_are_purely_real = all(imag(cross_corr_array(:)) == 0);

% Check whether correlation planes are complex.
if planes_are_purely_real
    % If planes are purely real, then automatically do the
    % displacement estimate in the spatial domain (peak finding)
   displacement_estimate_domain_string = 'spatial';
   ensemble_domain_string = 'spatial';
else
    % If planes are complex, then either spatial or
    % spectral displacement estimates will work.
    % In this case, determine the domain in which to
    % calculate displacement by reading the parameter
    % from the job file.
    displacement_estimate_domain_string = ...
        get_displacement_estimate_domain(JOBFILE, PASS_NUMBER);

    % Also set the "ensemble domain" specifier to "spectral"
    % so that the correct part of the switch/case statement
    % below gets chosen. 
    % This is a dumb way to handling this decision
    ensemble_domain_string = 'spectral';
end

% Spectral weighting method
spectral_weighting_method_string = ...
    lower(JOBFILE.Processing(PASS_NUMBER). ...
    Correlation.SpectralWeighting.Method);

% If the ensemble direction was "spatial,"
% then add all the planes together
switch lower(ensemble_direction_string)
    case 'spatial'
        % If the spatial correlation was specified,
        % then add all the correlations together. 
        cross_corr_array = sum(cross_corr_array, 3);
end

% Get the number of correlation planes remaining.
[region_height, region_width, num_correlation_planes] ...
    = size(cross_corr_array);

% Static particle diameters
% Particle diameter
particle_diameter = ...
   JOBFILE.Processing(PASS_NUMBER). ...
   Correlation.EstimatedParticleDiameter;

% Allocate the static spectral filter
% This apparently needs to be done 
% so that the parallel loop can run
spectral_weighting_filter_static = nan(region_height, region_width);

% Initialize a flag specifying whether
% spectral filters were pre-computed.
apc_do_temporal_ensemble = false;

% Filters were pre computed?
filters_were_pre_computed = false;

% Method specific options
switch lower(spectral_weighting_method_string)
    case 'rpc'   
        % Create the RPC filter
        spectral_weighting_filter_static = spectral_energy_filter( ...
            region_height, region_width, particle_diameter);
    case 'gcc'
    % Create the GCC filter (ones everywhere)
    spectral_weighting_filter_static = ones(region_height, region_width);
    
    case 'apc' 
        % If APC was selected, then check the APC method.
        % Read the APC method
        apc_method = JOBFILE.Processing(PASS_NUMBER). ...
            Correlation.SpectralWeighting.APC.Method;
        
        % Check whether the APC filters were pre computed.
        apc_do_temporal_ensemble = JOBFILE.Processing(PASS_NUMBER). ...
            Correlation.SpectralWeighting.APC.DoTemporalEnsemble;
        
        % Read them
        if apc_do_temporal_ensemble            
            % Read the pre-computed standard deviations
            apc_std_dev_x = ...
            	JOBFILE.Processing(PASS_NUMBER). ...
                Results.Filtering.APC.Diameter.X(:, 1);  
            apc_std_dev_y = ...
                JOBFILE.Processing(PASS_NUMBER). ...
                Results.Filtering.APC.Diameter.Y(:, 1);
            
            % Check whether any of the apc standard 
            % deviations were updated from NaN.
            % If not, then no filters were ever computed.
            filters_were_pre_computed = any(~isnan(apc_std_dev_x));
            
        end
end


% Make the list of particle diameters
% for the static methods
particle_diameter_list_x = zeros(num_correlation_planes, 1);
particle_diameter_list_y = zeros(num_correlation_planes, 1);

% Set the diameters
particle_diameter_list_x(:) = particle_diameter;
particle_diameter_list_y(:) = particle_diameter;

% Allocate arrays to hold vectors
TX = nan(num_correlation_planes, 1);
TY = nan(num_correlation_planes, 1);

% Allocate the subpixel weights
sub_pixel_weights = ones(region_height, region_width);

% Get SPC parameters if requested
switch lower(displacement_estimate_domain_string)
    case 'spectral'  
        % Get the list of kernel sizes for the SPC 
        spc_unwrap_method_string = ...
            get_spc_unwrap_method(JOBFILE, PASS_NUMBER);
        
        % Choose whether to run compiled codes for SPC
        spc_run_compiled = JOBFILE.Processing(PASS_NUMBER). ...
            Correlation.DisplacementEstimate. ...
            Spectral.RunCompiled;
        
end

% If the spectral ensemble was performed,
% then the spectal ensemble correlations
% need to be inverse-Fourier transformed
% back into the spatial domain. 
% This checks which type of ensemble was run
% (spatial or spectral) and does the inverse
% transform if necessary. Note that this transform
% happens in-place, i.e., the value of 
% the variable cross_corr_ensemble is changed
% after this result, rather than the transformed
% correlations being saved as a separate variable.
% This is to save memory, and also to reduce
% the number of variables we have to keep track of.
switch lower(ensemble_domain_string)
    case 'spectral'
        
        % Inform the user
        fprintf(1, 'Calculating inverse FTs...\n');
        % Do the inverse transform for each region.
        for k = 1 : num_correlation_planes
            
            % Extract the given region
            cross_corr_spectral = cross_corr_array(:, :, k);
            
            % Spectral correlation phase and magnitude
            [spectral_corr_phase, spectral_corr_mag] = ...
                split_complex(cross_corr_spectral);
                        
            % Switch between correlation methods
            switch lower(spectral_weighting_method_string)
                
                % SCC: Set the spectral magnitude as the weighting filter
                case 'scc'           
                    spectral_filter_temp = spectral_corr_mag;
                
                % APC: Use a dynamically-computed weighting filter.
                case 'apc'  
                    
                    % If the temporal ensemble was used to pre-compute
                    % the APC filter, then just read those values.
                    if apc_do_temporal_ensemble && filters_were_pre_computed
                        
                        % Read the horizontal standard deviation
                        % of the filter
                        filter_std_x = apc_std_dev_x(k);
                        filter_std_y = apc_std_dev_y(k);
                        
                        % Make the filter
                        spectral_filter_temp = ...
                            make_gaussian_spectral_filter(...
                            region_height, region_width, ...
                            filter_std_y, filter_std_x);
                    
                    % If the temporal ensemble wasn't used to pre-compute
                    % the APC filter, then compute it here.
                    else
                        % Calculate the APC filter
                        [spectral_filter_temp, filter_std_y, filter_std_x] = ...
                        calculate_apc_filter(cross_corr_spectral, ...
                        particle_diameter, apc_method);
                    end

                    % Equivalent particle diameter in the columns
                    % direction, calculated from the APC filter.
                    particle_diameter_list_x(k) = filter_std_dev_to_particle_diameter(...
                        filter_std_x, region_width);
                    
                    % Equivalent particle diameter in the rows
                    % direction, calculated from the APC filter.
                    particle_diameter_list_y(k) = filter_std_dev_to_particle_diameter(...
                        filter_std_y, region_height);
                otherwise
                    spectral_filter_temp = spectral_weighting_filter_static;
            end
            
            % Switch between estimating the 
            % displacement in the spatial domain (peak finding)
            % or in the spectral domain (plane fitting)
            switch lower(displacement_estimate_domain_string)
                case 'spatial'                    
                    % Apply the phase filter to the spectral plane.
                    % Note that this operation is legitimate even
                    % if SCC or GCC is selected. The reason for this
                    % is that if SCC is selected, then the "filter"
                    % is just the original magnitude. 
                    % I (Matt Giarra) have coded it this way to
                    % simplify the control flow. Not sure if this will
                    % turn out to be a nice way to do it.
                    cross_corr_filtered_spectral = ...
                    spectral_filter_temp .* spectral_corr_phase;

                    % Take the inverse FT of the "filtered" correlation.
                    cross_corr_filtered_spatial = fftshift(abs(ifft2(fftshift(...
                                cross_corr_filtered_spectral)))); 
                            
                    % Effective particle diameters
                    dp_x = particle_diameter_list_x(k);
                    dp_y = particle_diameter_list_y(k);

                    % Do the subpixel displacement estimate.
                    [TX(k), TY(k)] = ...
                        subpixel(cross_corr_filtered_spatial,...
                            region_width, region_height, sub_pixel_weights, ...
                                1, 0, [dp_x, dp_y]);              
                
                % If the displacement-estimate domain
                % was specified as "spectral" then
                % do the SPC plane fit to estimate the
                % pattern displacement.
                case 'spectral'    
                    % Calculate the displacement in the Fourier domain
                    % by unwrapping the phase angle and
                    % fitting a plane to it.
                    [TY(k), TX(k)] = spc_2D(...
                        cross_corr_spectral, ...
                        spectral_filter_temp, ...
                        spc_unwrap_method_string, ...
                        spc_run_compiled);
            end                 
        end
        
    % This is the case for ensembles done in the spatial
    % domain rather than the spectral domain. This case
    % assumes that the correlation planes passed are
    % purely real. 
    otherwise
        
        % If the displacement estimate was specified
        % to occur in the spatial domain (peak finding)
        % then do the subpixel peak fitting. 
        % For now the three-point Gaussian fit is hard coded.
        % This should be changed.
        switch lower(displacement_estimate_domain_string)
            case 'spatial'
                % Loop over the planes
                for k = 1 : num_correlation_planes
                    cross_corr_spatial = cross_corr_array(:, :, k);

                    % Effective particle diameters
                    dp_x = particle_diameter_list_x(k);
                    dp_y = particle_diameter_list_y(k);

                    % Do the subpixel displacement estimate.
                    [TX(k), TY(k)] = ...
                        subpixel(cross_corr_spatial,...
                        region_width, region_height, sub_pixel_weights, ...
                        1, 0, [dp_x, dp_y]);   
                end
        end
end

% Save the particle diameters
PARTICLE_DIAMETER_Y = particle_diameter_list_y;
PARTICLE_DIAMETER_X = particle_diameter_list_x;

end


