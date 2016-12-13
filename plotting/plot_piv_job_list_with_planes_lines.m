function plot_piv_job_list_with_planes_lines(JOBLIST, PASS_NUMBER, VECT_TO_PLOT)

% Add paths
addpaths();

% Font size
fSize_textbox = 20;

% Font size for the colorbar
fSize_colorbar = 20;

fSize_title = 12;

% Vector scale
% Dense plot:
vect_scale = 4;

% Vector line width
% Sparse plot: 
% vector_line_width = 0.1;

% Dense plot:
vector_line_width = 0.5;

% Vector skip

% Sparse plot: 
skip_x = 10;
skip_y = 2;


% Dense plot:
% skip_x = 4;
% skip_y = 4;

% Vector offset
x_offset = 1;
y_offset = 0;

% Plane X and Y spots
gx_plane = 1792;
gy_plane = 657;

gx_plane = 1664;
gy_plane = 689;

% gx_plane = 1680;
% gy_lane = 689;

plane_skip_x = 3;
plane_skip_y = 3;

yc = 580;
xc = 1600;

plot_width = 550;

xl_sub = xc + (round(plot_width/2) * [-1, 1]);
yl_sub = yc + (round(plot_width/2) * [-1, 1]);

% Color scale
color_scale = [-20, 50];

% Subplot margins
subplot_margin_width  = 0.1 * [1, 1.5];
subplot_margin_height = 0.1 * [0.1 1];
% subplot_gap = 0.01 * [0, 1];
subplot_gap = [];

% Plot vertical shift
dy_subplot = 60;

% Plot horizontal shift;
dx_subplot = -45;

dx_meshplot = -30;

% Location of the bottom edge of the top-left plot
ax_pos_bottom_main = 245;

% Height of the top left plot
ax_height_main = 177;

% Default to plotting whatever
% post processing was done
% rather than a specific
% step (e.g., raw, validated, smoothed)
if nargin < 3
    VECT_TO_PLOT = 'final';
end

% Default to plotting the final pass
if nargin < 2
    PASS_NUMBER = 0;
end


% Set pass number to 0 if an empty argument is passed.
if isempty(PASS_NUMBER)
    PASS_NUMBER = 0;
end

plot_final = ~isempty(regexpi(VECT_TO_PLOT, 'fi'));
plot_validated = ~isempty(regexpi(VECT_TO_PLOT, 'va'));
plot_smoothed = ~isempty(regexpi(VECT_TO_PLOT, 'sm'));
plot_raw = ~isempty(regexpi(VECT_TO_PLOT, 'raw'));


% % % % % % % % %


% Count the number of jobs
num_jobs = length(JOBLIST);

% Initialize the max and min velocities.
vel_max = -inf;
vel_min = inf;

% Axes' positions
ax_pos = zeros(num_jobs, 4);

% Make a figure
figure(1);

% Delete annotations
delete(findall(gcf,'Tag','label_tag'))

% Number of correlations
corr_strings = {};

% Image base names
image_base_names = {};

% Determine which correlations
% and images were correlated.
% Caveat: This assumes that
% the correlation type was 
% the same for all passes.
for n = 1 : num_jobs
    
    JobFile = JOBLIST(n);
    
    % Correlation type
    corr_strings{n} = JobFile.Processing(1).Correlation.Method;

    % Number of types of images
    image_base_names{n} = JobFile.Data.Inputs.Images.BaseName;  
    
end

% Number of different correlation methods
num_correlation_methods = length(unique(corr_strings));

% Subplot rows and cols
subplot_rows = num_jobs;
subplot_cols = 3;

% Loop over the jobs
for n = 1 : num_jobs
    
    % Extract the job file
    JobFile = JOBLIST(n);
    
    % Determine the save path
    output_file_path = determine_jobfile_save_path(JobFile);

    % Sub job with just the little
    % region done and containing planes.
    sub_job_path = strrep(output_file_path, 'c_0', 'c_test_0');
    
    % Load the output file
    loaded_jobfile = load(output_file_path);
    
    % Loaded subjob
    loaded_subjobfile = load(sub_job_path);
    
    % Loaded job file
    JobFile_Loaded = loaded_jobfile.JobFile;
    
    % Measure the number of passes
    num_passes = length(JobFile_Loaded.Processing);
    
    % Figure out which pass to plot
    if PASS_NUMBER == 0
        pass_to_plot = num_passes;
    else
        % Figure out which pass to plot
        pass_to_plot = min(PASS_NUMBER, num_passes);
    end
    
    % Determine which displacement field to plot
    if plot_raw
        displacement_field = JobFile_Loaded.Processing(pass_to_plot).Results.Displacement.Raw;
    elseif plot_validated
        displacement_field = JobFile_Loaded.Processing(pass_to_plot).Results.Displacement.Validated;
    elseif plot_smoothed
        displacement_field = JobFile_Loaded.Processing(pass_to_plot).Results.Displacement.Smoothed;
    elseif plot_final
        displacement_field = JobFile_Loaded.Processing(pass_to_plot).Results.Displacement.Final;
    end

    % Read the displacement field
    tx = displacement_field.X(:, 1);
    ty = displacement_field.Y(:, 1);
    
    % Read the grid points.
    gx = JobFile_Loaded.Processing(pass_to_plot).Grid.Points.Full.X;
    gy = JobFile_Loaded.Processing(pass_to_plot).Grid.Points.Full.Y;
    
    % Count the grid points
    nx = length(unique(gx(:)));
    ny = length(unique(gy(:)));
    
    % Reshape grid
    gx_mat = reshape(gx, [ny, nx]);
    gy_mat = reshape(gy, [ny, nx]);
    
    % Reshape vectors
    tx_mat = reshape(tx, [ny, nx]);
    ty_mat = reshape(ty, [ny, nx]);
    
    % Subplot row
    subplot_row = n;
    
    % Create a figure
    subtightplot(subplot_rows, subplot_cols, ...
        subplot_cols * (n - 1) + 1, ...
        subplot_gap, ...
        subplot_margin_height, ...
        subplot_margin_width);
    set(gca, 'units', 'pixels');
   
    % Make a plot
    % This will probably screw up the positions?
    imagesc(gx, gy, tx_mat);
    axis image;
    hold on;
    quiver(             gx_mat((1 + y_offset) : skip_y : end, (1 + x_offset) : skip_x : end), ...
                        gy_mat((1 + y_offset) : skip_y : end, (1 + x_offset) : skip_x : end), ...
           vect_scale * tx_mat((1 + y_offset) : skip_y : end, (1 + x_offset) : skip_x : end), ...
           vect_scale * ty_mat((1 + y_offset) : skip_y : end, (1 + x_offset) : skip_x : end), ...
           0, '-k', 'linewidth', vector_line_width);
       
    % Format axes
    axis off
   
    % Load the mask
    grid_mask = load_mask(JobFile, num_passes);
    [mask_points_y, mask_points_x] = get_mask_outline(grid_mask);

    % Plot the mask points.
    plot(mask_points_x, mask_points_y, 'ow', 'markersize', 1, 'markerfacecolor', 'white');
    
    % Show the subregion.
    % rect_coordinates
    rect_y = [yl_sub(1), yl_sub(1), yl_sub(2), yl_sub(2), yl_sub(1)];
    rect_x = [xl_sub(1), xl_sub(2), xl_sub(2), xl_sub(1), xl_sub(1)];
    h_rect = plot(rect_x, rect_y, '--w');
    
    plot(gx_plane, gy_plane, 'ow', 'markerfacecolor', 'w', 'markersize', 1.5);
    
    
    hold off;    
    
    % Calculate the max and min velocities
    % of all the plotted data.
    vel_max = max(vel_max, max(tx(:)));
    vel_min = min(vel_min, min(tx(:)));

    % Set the color scale
    if isempty(color_scale)
        caxis([vel_min, vel_max]);
    else
        caxis(color_scale);
    end
 
    % Plot limits
    xlim([min(gx(:)), max(gx(:))]);
    ylim([min(gy(:)), max(gy(:))]);

   
    % Move the plot vertically 
    % because subtightplot() seems
    % to not be able to make the vertical
    % gap small enough for my preferences.
    ax_pos_main_current = get(gca, 'outerposition');
    ax_pos_main_bottom = ax_pos_main_current(2) + (n - 1) * dy_subplot;
    ax_pos_main_current(2) = ax_pos_main_bottom;
%     ax_pos_main_current(2) = ax_bottom_main + ...
%         (n - 1) * (dy_subplot - ax_height_main);
    ax_pos_main_current(4) = ax_height_main;
   
    % Position of the right edge of this plot
    ax_right_main_current(1) = ax_pos_main_current(1) + ...
        ax_pos_main_current(3); 
    
     set(gca, 'outerposition', ax_pos_main_current);
    
    
%     % Location of bottom of axes.
%     ax_bottom_main = ax_pos_main_current(2);
%     
%     % Height of axes.
%     ax_height_main = ax_pos_main_current(4);
    
    
    % Flip the axis direction
    set(gca, 'ydir', 'normal');  

    % Axes' positions
    ax_pos(n, :) = get(gca, 'position');
    

   
    
    
    % Second subplot % % % %
    % % % % % % % %
    
    gx_val = gx(gx >= xl_sub(1) & gx <= xl_sub(2) & ...
                gy >= yl_sub(1) & gy <= yl_sub(2));

    gy_val = gy(gx >= xl_sub(1) & gx <= xl_sub(2) & ...
        gy >= yl_sub(1) & gy <= yl_sub(2));
    
    tx_val = tx(gx >= xl_sub(1) & gx <= xl_sub(2) & ...
                gy >= yl_sub(1) & gy <= yl_sub(2));

    ty_val = ty(gx >= xl_sub(1) & gx <= xl_sub(2) & ...
        gy >= yl_sub(1) & gy <= yl_sub(2));
    
    % Number in each direction
    nx_val = length(unique(gx_val));
    ny_val = length(unique(gy_val));
    
    % Reshape displacements into grids
    tx_val_mat = reshape(tx_val, [ny_val, nx_val]);
    ty_val_mat = reshape(ty_val, [ny_val, nx_val]);
    
    % average along X
    tx_val_mean = mean(tx_val_mat, 2);
    ty_val_mean = mean(ty_val_mat, 2);
    
    tx_val_std = std(tx_val_mat, [], 2);
    ty_val_std = std(ty_val_mat, [], 2);

    % Create a subplot
    subtightplot(subplot_rows, subplot_cols, ...
        subplot_cols * (n - 1) + 2, ...
        subplot_gap, ...
        subplot_margin_height, ...
        subplot_margin_width);
    set(gca, 'units', 'pixels');
    
    % Unique Y grid points
    y_subplot = unique(gy_val);
    
    % X limits   
    xlim_subplot = [y_subplot(1), y_subplot(end)];

    shadedErrorBar(unique(gy_val), tx_val_mean, tx_val_std);
%     shadedErrorBar([], tx_val_mean, tx_val_std)
    axis square;
    set(gca, 'xdir', 'reverse');
    set(gca, 'view', [90, 90]);
    xlim(xlim_subplot);
    ylim([-18, 45]);
    set(gca,'xaxisLocation','top')
    
    % Move the plot vertically 
    % because subtightplot() seems
    % to not be able to make the vertical
    % gap small enough for my preferences.
    ax_pos_subplot_current = get(gca, 'outerposition');

    % Left edge
    ax_pos_subplot_current(1) = ax_pos_subplot_current(1) + 1 * dx_subplot;
    
    % Bottom edge
    ax_pos_subplot_current(2) = ax_pos_main_bottom;
    
    % Height (dont need to set width because
    % aspect ratio is already fixed)
    ax_pos_subplot_current(4) = ax_height_main;
    
    % Set the plot position
    set(gca, 'outerposition', ax_pos_subplot_current);
    
    % Flip the axis direction
    set(gca, 'ydir', 'normal');
    
    if n < num_jobs
        set(gca, 'yticklabel', '');
    end
    
    if n == 1
        title('$\Delta x \left( \textrm{pixels} \right)$', ...
            'interpreter', 'latex', 'fontsize', fSize_title);
    end
    
    % Turn off the other ticks
    set(gca, 'xticklabel', '');
    
    
    % % % % % % % % % % % %
    % % % % Third Subplot
    % % % % % % % % % % % %
    % Create a subplot
    subtightplot(subplot_rows, subplot_cols, ...
        subplot_cols * (n - 1) + 3, ...
        subplot_gap, ...
        subplot_margin_height, ...
        subplot_margin_width);
    set(gca, 'units', 'pixels');
    
    % Extract the planes
    c_planes   = loaded_subjobfile.JobFile.Processing(1).Results.Planes;
    gx_subplot = loaded_subjobfile.JobFile.Processing(1).Grid.Points.Correlate.X;
    gy_subplot = loaded_subjobfile.JobFile.Processing(1).Grid.Points.Correlate.Y;
    
    % Region size
    [region_height, region_width, num_planes] = size(c_planes);
    
    % Find the index corresponding to
    % the specified coordinates of the plane.
    plane_ind = find(gx_subplot == gx_plane & gy_subplot == gy_plane);

    % Extract the plane
    cplane = c_planes(1 : plane_skip_y : end, 1 : plane_skip_x : end, plane_ind);
    cplane_sub = cplane - min(cplane(:));
    
    cplane_norm = cplane_sub ./ max(cplane_sub(:));
    
    % Mesh plot
    mesh(cplane_norm, 'edgecolor', 'black', 'linewidth', 0.1);
    axis square;
    xlim([1, region_width/plane_skip_x]);
    ylim([1, region_height / plane_skip_y]);
    zlim([0, 1])
    set(gca, 'view', [-20, 10]);
   
    p = get(gca, 'outerposition');
    p(1) = p(1) + 2 * dx_subplot + dx_meshplot;
    p(2) = ax_pos_main_bottom;
    p(4) = ax_height_main;
    set(gca, 'outerposition', p);
    axis off;
    
      
end

set(gcf, 'color' ,'white');

fprintf(1, '');



% % % % 
% % % % % % Make a color bar
% % % % h_colorbar = colorbar;
% % % % set(h_colorbar, 'units', 'pixels');
% % % % 
% % % % % This is the axis position
% % % % % of the last axis that was plotted.
% % % % ax_pos_end = ax_pos(end, :);
% % % % 
% % % % % Reset the axis position
% % % % set(gca, 'position', ax_pos_end);
% % % % 
% % % % dx_cbar = 3;
% % % % 
% % % % % Full plots: use 57
% % % % % Subplots: use ....
% % % % dy_cbar = 57;
% % % % 
% % % % % Full plots: Use 0.817
% % % % height_fract_cbar = 0.817;
% % % % 
% % % % % Figure out the location of
% % % % % the bottom edge of the bottom-most plot
% % % % plot_bottom_edge = ax_pos(end, 2);
% % % % 
% % % % % Figure out the location of
% % % % % the top edge of the top-most subplot
% % % % plot_top_edge = ax_pos(1, 2) + ax_pos(1, 4);
% % % % 
% % % % % Plot right edge
% % % % plot_right_edge = ax_pos(end, 1) + ax_pos(end, 3);
% % % % 
% % % % % colorbar_position
% % % % pos_colorbar = get(h_colorbar, 'position');
% % % % 
% % % % % Shift the color bar horizontally
% % % % pos_colorbar(1) = plot_right_edge + dx_cbar;
% % % % 
% % % % % Set the bottom of the color bar
% % % % % to be aligned with the bottom
% % % % % edge of the bottom plot
% % % % pos_colorbar(2) = plot_bottom_edge + dy_cbar;
% % % % 
% % % % % Set the height of the color
% % % % % bar to be equal to the distance
% % % % % between the top edge of the top plot
% % % % % and the bottom edge of the bottom plot.
% % % % pos_colorbar(4) = height_fract_cbar * (plot_top_edge - plot_bottom_edge);
% % % % 
% % % % % Update the position of the color bar
% % % % set(h_colorbar, 'position', pos_colorbar);
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % % Set colorbar tick font size
% % % % set(h_colorbar, 'fontsize', 0.8 * fSize_colorbar);
% % % % 
% % % % % Set color bar label
% % % % ylabel(h_colorbar, '$\Delta X \, \left( \textrm{pixels} \right)$', ...
% % % %     'interpreter', 'latex', 'fontsize', fSize_colorbar);
% % % % 
% % % % 
% % % % 
% % % % % % % 
% % % % %   %
% % % % %   %
% % % % % % % 
% % % % 
% % % % 
% This block of code creates the labels on the 
% plots that indicate the types of
% correlations that were performed. 
% Unique correlation strings

% Delete any existing annotations.
delete(findall(gcf,'Tag','label_tag'))

% Extract the unique strings for correlation type.
% unique_corr_strings = fliplr((unique(corr_strings)));


% These two variables control
% the positions of the labels
dx_label_corr_type = -45;
dy_label_corr_type = -5;
for n = 1 : num_correlation_methods
    
    % Position vector
    pos_vect = ax_pos(n, :);
    
    % Vertical shift of the label
    dy_label = dy_label_corr_type + pos_vect(4)/2;
    
    % Label position
    label_pos_x = pos_vect(1) + dx_label_corr_type;
    label_pos_y = pos_vect(2) + dy_label;
    
    % Label position
    label_pos = [label_pos_x, label_pos_y, 0, 0];
    
    % Make the annotation
    h = annotation('textbox', ...
    'string', upper(corr_strings{n}), ...
    'fontsize', fSize_textbox, ...
    'interpreter', 'latex', ...
    'Tag', 'label_tag', ...
    'linestyle', 'none',...
    'units', 'pixels', ...
    'verticalalignment', 'middle', ...
    'position', label_pos);    
end
% % % % 
% % % % % % % 
% % % % %   %
% % % % %   %
% % % % % % % 
% % % % 
% % % % 
% % % % 
% % % % % % % 
% % % % %   %
% % % % %   %
% % % % % % % 
% % % % 
% % % % 
% % % % % This block of code creates the labels on the 
% % % % % plots that indicate the types of
% % % % % images were correlated. 
% % % % 
% % % % % Extract the unique strings
% % % % % for the image base names.
% % % % unique_image_names = unique(image_base_names);
% % % % 
% % % % % These two variables control
% % % % % the positions of the labels
% % % % dx_label_image_names = 0;
% % % % dy_label_image_names = 20;
% % % % for n = 1 : num_image_types
% % % %     % Axis index
% % % %     ax_ind = sub2ind([subplot_rows, subplot_cols], 1, n);
% % % %     
% % % %     % Label string
% % % %     img_type_current = unique_image_names{n};
% % % %     
% % % %     % Determine what string to plot.
% % % %     % We don't want to plot the actual
% % % %     % image names since that makes
% % % %     % for a really dumb looking plot.
% % % %     if isempty(regexpi(img_type_current, 'gho'))
% % % %         label_string = 'Raw';
% % % %     else
% % % %         label_string = 'Pre-processed';
% % % %     end
% % % %     
% % % %     % Latex formatted label string
% % % %     label_string_latex = sprintf('\\textrm{%s}', label_string);
% % % %     
% % % %     % Position vector for the top row axis.
% % % %     pos_vect = ax_pos(ax_ind, :);
% % % %     
% % % %     % Axis x and y positions
% % % %     ax_pos_left = pos_vect(1);
% % % %     ax_pos_bottom = pos_vect(2);
% % % %     ax_width = pos_vect(3);
% % % %     ax_height = pos_vect(4);
% % % %     
% % % %     % Vertical shift of the label
% % % %     dy_label = dy_label_image_names + ax_height;
% % % %     dx_label = dx_label_image_names + ax_width  / 2;
% % % %     
% % % %     label_pos_x = ax_pos_left   + dx_label;
% % % %     label_pos_y = ax_pos_bottom + dy_label;
% % % %     
% % % %     % Label position
% % % %     label_pos = [label_pos_x, label_pos_y, 0, 0];
% % % %     
% % % %     % Make the annotation
% % % %     h = annotation('textbox', ...
% % % %     'string', label_string_latex, ...
% % % %     'fontsize', fSize_textbox, ...
% % % %     'interpreter', 'latex', ...
% % % %     'Tag', 'label_tag', ...
% % % %     'linestyle', 'none',...
% % % %     'units', 'pixels',...
% % % %     'horizontalAlignment', 'center', ...
% % % %     'position', label_pos);
% % % % 
% % % % %     set(h, 'position', label_pos);
% % % %     
% % % % end
% % % % 
% % % % % % % 
% % % % %   %
% % % % %   %
% % % % % % % 
% % % % 
% % % % % Make the figure white
% % % % set(gcf, 'color', 'white');
% % % % 
% % % % % % Make the labels
% % % % % for n = 1 : num_jobs
% % % % %    
% % % % %     if corr_strings{n} = unique_corr
% % % % %     
% % % % % end
% % % % 
% % % % 
% % % % % % Textbox position
% % % % % textbox_pos = ax_pos(1, :);
% % % % % 
% % % % % % textbox_x = 
% % % % % 
% % % % % % Annotations
% % % % % annotation('textbox', 'position', textbox_pos, ...
% % % % %     'string', 'Fit', ...
% % % % %     'fontsize', fSize_textbox, ...
% % % % %     'interpreter', 'latex', ...
% % % % %     'Tag', 'legend_text_01', 'linestyle', 'none');


end





