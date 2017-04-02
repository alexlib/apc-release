function GAUSSIAN_SPECTRAL_FILTER = make_gaussian_spectral_filter(...
    REGION_HEIGHT, REGION_WIDTH, GAUSSIAN_STD_DEV_Y, GAUSSIAN_STD_DEV_X)

    % Default to symmetric filter
    if nargin < 4
        GAUSSIAN_STD_DEV_X = GAUSSIAN_STD_DEV_Y;
    end

    % Make coordinate vectors
    xv = (1 : REGION_WIDTH) - fourier_zero(REGION_WIDTH);
    yv = (1 : REGION_HEIGHT) - fourier_zero(REGION_HEIGHT);
    
    % Make coordinate arrays
    [X, Y] = meshgrid(xv, yv);
    
    % Calculate filter
    GAUSSIAN_SPECTRAL_FILTER = exp(-X.^2 / (2 * GAUSSIAN_STD_DEV_X^2) - Y.^2 / (2 * GAUSSIAN_STD_DEV_Y^2));
    
end