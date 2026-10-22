function plotChannel(r_mat, fluid, options)
% PLOTCHANNEL Generates docked figures analyzing a rectangular cooling channel.
%
% Inputs:
%   r_mat   - [2 x N+1] Nodal coordinates [x; y] of the channel path [m].
%   fluid   - [struct] Contains thermophysical and geometric properties evaluated at midpoints.
%   options - [name-value] Optional arguments (e.g., title_prefix).

arguments
    r_mat (2,:) double
    fluid struct
    options.title_prefix (1,1) string = ""
end

    prefix = "";
    if options.title_prefix ~= ""
        prefix = options.title_prefix + " - ";
    end

    num_nodes = size(r_mat, 2);
    x_axis = r_mat(1, 2:end);
    
    fig_name_3d = prefix + "3D Extruded Fluid Channel";
    fig3d = figure('Name', fig_name_3d, 'WindowStyle', 'docked');
    ax3d = axes(fig3d);
    hold(ax3d, 'on');
    
    dx = diff(r_mat(1, :));
    dy = diff(r_mat(2, :));
    ds = sqrt(dx.^2 + dy.^2);
    s_nodes = [0, cumsum(ds)];          
    s_mid = s_nodes(1:end-1) + ds / 2;  
    
    T_nodes = interp1(s_mid, fluid.Tprofile, s_nodes, 'linear', 'extrap');
    a_nodes = interp1(s_mid, fluid.ChannelBase, s_nodes, 'linear', 'extrap');
    b_nodes = interp1(s_mid, fluid.ChannelHeight, s_nodes, 'linear', 'extrap');
    
    r_3d = [r_mat; zeros(1, num_nodes)]; 
    T_vec = zeros(3, num_nodes);
    for i = 1:num_nodes
        if i == 1
            dir = r_3d(:, 2) - r_3d(:, 1);
        elseif i == num_nodes
            dir = r_3d(:, end) - r_3d(:, end-1);
        else
            dir = r_3d(:, i+1) - r_3d(:, i-1); 
        end
        T_vec(:, i) = dir / norm(dir);
    end
    
    N_vec = [-T_vec(2, :); T_vec(1, :); zeros(1, num_nodes)];
    B_vec = repmat([0; 0; 1], 1, num_nodes);
    
    num_points = 5;
    X = zeros(num_nodes, num_points);
    Y = zeros(num_nodes, num_points);
    Z = zeros(num_nodes, num_points);
    C = zeros(num_nodes, num_points);
    
    for i = 1:num_nodes
        a = a_nodes(i);
        b = b_nodes(i);
        
        corners_x = [ (b/2)*N_vec(1,i) + (a/2)*B_vec(1,i), ...
                      (b/2)*N_vec(1,i) - (a/2)*B_vec(1,i), ...
                     -(b/2)*N_vec(1,i) - (a/2)*B_vec(1,i), ...
                     -(b/2)*N_vec(1,i) + (a/2)*B_vec(1,i), ...
                      (b/2)*N_vec(1,i) + (a/2)*B_vec(1,i) ]; 
                  
        corners_y = [ (b/2)*N_vec(2,i) + (a/2)*B_vec(2,i), ...
                      (b/2)*N_vec(2,i) - (a/2)*B_vec(2,i), ...
                     -(b/2)*N_vec(2,i) - (a/2)*B_vec(2,i), ...
                     -(b/2)*N_vec(2,i) + (a/2)*B_vec(2,i), ...
                      (b/2)*N_vec(2,i) + (a/2)*B_vec(2,i) ];
                  
        corners_z = [ (b/2)*N_vec(3,i) + (a/2)*B_vec(3,i), ...
                      (b/2)*N_vec(3,i) - (a/2)*B_vec(3,i), ...
                     -(b/2)*N_vec(3,i) - (a/2)*B_vec(3,i), ...
                     -(b/2)*N_vec(3,i) + (a/2)*B_vec(3,i), ...
                      (b/2)*N_vec(3,i) + (a/2)*B_vec(3,i) ];
        
        X(i, :) = r_3d(1, i) + corners_x;
        Y(i, :) = r_3d(2, i) + corners_y;
        Z(i, :) = r_3d(3, i) + corners_z;
        C(i, :) = T_nodes(i); 
    end
    
    surf(ax3d, X, Y, Z, C, 'EdgeColor', '#444444', 'FaceColor', 'interp');
    colormap(ax3d, 'turbo');
    cb = colorbar(ax3d);
    cb.Label.String = 'Bulk Fluid Temperature [K]';
    
    axis(ax3d, 'equal'); 
    grid(ax3d, 'on'); 
    view(ax3d, 3);
    camlight(ax3d, 'headlight'); 
    lighting(ax3d, 'gouraud'); 
    material(ax3d, 'shiny');
    
    xlabel(ax3d, 'Axial Position, x [m]');
    ylabel(ax3d, 'Radial Position, y [m]');
    zlabel(ax3d, 'Circumferential, z [m]');
    title(ax3d, prefix + "Continuous Rectangular Channel", 'FontWeight', 'bold');
    
    generate_docked_plot(x_axis, fluid.Tprofile,          prefix, 'Bulk Fluid Temperature', 'T [K]', '#0677c2');
    generate_docked_plot(x_axis, fluid.HcoefConv,         prefix, 'Heat Transfer Coefficient', 'h [W/m^2K]', '#d95319');
    generate_docked_plot(x_axis, fluid.Velocity,          prefix, 'Fluid Velocity', 'v [m/s]', '#edb120');
    
    generate_docked_plot(x_axis, fluid.ChannelBase,       prefix, 'Channel Base Width (a)', 'a [m]', '#7e2f8e');
    generate_docked_plot(x_axis, fluid.ChannelHeight,     prefix, 'Channel Height (b)', 'b [m]', '#7e2f8e');
    generate_docked_plot(x_axis, fluid.HydraulicDiameter, prefix, 'Hydraulic Diameter (Dh)', 'D_h [m]', '#7e2f8e');
    generate_docked_plot(x_axis, fluid.BaseOverHeight,    prefix, 'Channel Aspect Ratio (a/b)', 'Aspect Ratio [-]', '#7e2f8e');
    
    generate_docked_plot(x_axis, fluid.Re,                prefix, 'Reynolds Number', 'Re [-]', '#77ac30');
    generate_docked_plot(x_axis, fluid.Pr,                prefix, 'Prandtl Number', 'Pr [-]', '#4dbeee');
    generate_docked_plot(x_axis, fluid.cp_vec,            prefix, 'Specific Heat Capacity', 'cp [J/kgK]', '#a2142f');
    generate_docked_plot(x_axis, fluid.mu_vec,            prefix, 'Dynamic Viscosity', 'mu [Pa*s]', '#0072bd');
    generate_docked_plot(x_axis, fluid.K_vec,             prefix, 'Thermal Conductivity', 'k [W/mK]', '#d95319');
    
    generate_docked_plot(x_axis, fluid.massFlowRateError, prefix, 'Mass Flow Rate Error', '\Delta mdot [kg/s]', '#000000');

end

function generate_docked_plot(x_data, y_data, prefix, prop_name, y_label_str, line_color)
% GENERATE_DOCKED_PLOT Helper function to standardize 2D plotting formatting.
%
% Inputs:
%   x_data      - [1 x N] Independent variable vector (axial position) [m].
%   y_data      - [1 x N] Dependent variable vector to be plotted.
%   prefix      - [string] Title prefix for the figure.
%   prop_name   - [string] Name of the property for the title.
%   y_label_str - [string] Y-axis label including dimensions.
%   line_color  - [char/string] Hexadecimal color code or standard MATLAB color character.

    full_title = prefix + prop_name;
    fig = figure('Name', full_title, 'WindowStyle', 'docked');
    ax = axes(fig);
    
    if startsWith(line_color, '#') && length(line_color) == 2
        line_color = line_color(2);
    end
    
    hLine = plot(ax, x_data, y_data, '-o', 'LineWidth', 1.5, 'Color', line_color, 'MarkerSize', 5);
    hLine.MarkerFaceColor = hLine.Color; 
    hLine.MarkerEdgeColor = 'none';
    
    title(ax, full_title, 'FontWeight', 'bold');
    xlabel(ax, 'Axial Position, x [m]');
    ylabel(ax, y_label_str, 'Interpreter', 'tex'); 
    grid(ax, 'on');
    set(ax, 'FontSize', 10);
end