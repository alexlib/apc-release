function [APC_FILTER, APC_STD_Y, APC_STD_X] = ...
    calculate_apc_filter_magnitude_method(SPECTRAL_CORRELATION_COMPLEX, RPC_DIAMETER)

    % Default to no minimum rpc diameter
    if nargin < 2
        RPC_DIAMETER = 0;
    end

    % Size of the region
    [region_height, region_width] = size(SPECTRAL_CORRELATION_COMPLEX);

    % Standard deviations of the RPC filter.
    rpc_std_dev_x = particle_diameter_to_filter_std_dev(RPC_DIAMETER, region_width);
    rpc_std_dev_y = particle_diameter_to_filter_std_dev(RPC_DIAMETER, region_height);
    
    % Fit a Gaussian function to the magnitude
    % of the complex correlation, 
    % which should represent the SNR versus wavenumber.
    [~, sy, sx] =...
        fit_gaussian_2D(abs(SPECTRAL_CORRELATION_COMPLEX) ...
        ./ max(abs(SPECTRAL_CORRELATION_COMPLEX(:))));
    
    % The fit can crap out and come back with
    % a standard deviation of less than 1. This is nonphysical
    % and can be used as a flag.
    if sx <= 1
        sx = rpc_std_dev_x;
    end
    if sy <= 1
        sy = rpc_std_dev_y;
    end
    
    % Take the APC diameter as the minimum
    % between the RPC equivalent std dev
    % and the standard deviation diameter.
    APC_STD_Y = min(rpc_std_dev_y, sy);
    APC_STD_X = min(rpc_std_dev_x, sx);
    
    % Calculate the filter
    APC_FILTER = make_gaussian_spectral_filter(...
        region_height, region_width, APC_STD_Y, APC_STD_X);
    
end





