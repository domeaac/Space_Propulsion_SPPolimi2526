function plotSolution(Solution, options)
% PLOTSOLUTION Generates axial and spatial temperature contour plots,
% aerodynamic, and aerothermal profiles in docked figures.
%
% Inputs:
%   Solution - Struct output from the steady-state thermal solver.
%   options.thicknessScales - (Optional) Scalar or vector specifying the 
%                             magnification factor for layer thicknesses. 
%                             If a vector, length must equal nLayers.
%   options.sectionTitle    - (Optional) String or char scalar to prepend 
%                             to the title of all docked figure tabs.

    arguments
        Solution (1,1) struct
        options.thicknessScales double = []
        options.sectionTitle {mustBeTextScalar} = ""
    end

    Mesh = Solution.Mesh;
    nSegments = Mesh.nCurveNodes - 1;
    nLayers = Mesh.nLayers;
    nNodes = Mesh.nCurveNodes;

    if isempty(options.thicknessScales)
        thicknessScales = ones(1, nLayers);
    elseif isscalar(options.thicknessScales)
        thicknessScales = options.thicknessScales * ones(1, nLayers);
    else
        assert(length(options.thicknessScales) == nLayers, ...
            'Length of thicknessScales vector must match nLayers.');
        thicknessScales = options.thicknessScales;
    end
    
    if options.sectionTitle ~= ""
        prefix = char(options.sectionTitle + " - ");
    else
        prefix = '';
    end
    
    x_mids = Solution.MidPointCoords(:, 1)';
    T_matrix = flipud(Solution.TemperatureProfile(:, 2:end)');
    
    T_aw = Solution.TemperatureProfile(:, 1)';
    q_flux = Solution.HeatFluxProfile'; 
    h_conv = Solution.ConvectionCoef';
    M_num  = Solution.Mach'; 
    
    TotalTemp = Solution.TotalTemperature(:)';
    cp        = Solution.cp(:)';
    gamma     = Solution.gamma(:)';
    prandtl   = Solution.prandtl(:)';
    mu        = Solution.mu(:)';
    T_cea     = Solution.T_cea(:)';
    AxAt      = Solution.AxAt(:)';
    cstar     = Solution.cstar(:)';
    
    if isscalar(cstar)
        cstar = cstar * ones(size(x_mids));
    end

    X_int = zeros(nNodes, nLayers + 1);
    Y_int = zeros(nNodes, nLayers + 1);
    
    c_nodes = Mesh.curveNodes;
    l_nodes = Mesh.layerNodes;
    
    X_int(:, 1) = l_nodes(1, :, 1)';
    Y_int(:, 1) = l_nodes(2, :, 1)';

    for k = 1:nLayers
        dx_layer = l_nodes(1, :, k+1) - l_nodes(1, :, k);
        dy_layer = l_nodes(2, :, k+1) - l_nodes(2, :, k);
        
        X_int(:, k+1) = X_int(:, k) + (thicknessScales(k) * dx_layer)';
        Y_int(:, k+1) = Y_int(:, k) + (thicknessScales(k) * dy_layer)';
    end

    figure('Name', [prefix 'Axial Temperature Distribution'], 'Color', 'w', 'WindowStyle', 'docked');
    hold on; grid on; box on;
    
    plot(x_mids, T_aw, 'k--', 'LineWidth', 1.5, 'DisplayName', 'T_{aw} (Adiabatic Wall)');
    colors = lines(nLayers + 1);
    
    for k = (nLayers + 1):-1:1
        if k == nLayers + 1
            displayName = 'Hot Inner Wall';
        elseif k == 1
            displayName = 'Cold Outer Wall';
        else
            displayName = sprintf('Layer Interface %d', k);
        end
        
        p = plot(x_mids, T_matrix(k, :), '-o', 'LineWidth', 1.5, 'Color', colors(k, :), ...
                 'MarkerSize', 4, 'DisplayName', displayName);
        p.MarkerFaceColor = p.Color;
    end
    
    xlabel('Axial Position, x (m)');
    ylabel('Temperature, T (K)');
    title('Axial Temperature Profile per Layer Interface');
    legend('Location', 'best');
    
    figure('Name', [prefix '2D Thermal Colorgrading'], 'Color', 'w', 'WindowStyle', 'docked');
    ax2 = gca;
    hold(ax2, 'on'); 
    axis(ax2, 'equal'); 
    box(ax2, 'on');
    
    ax2.SortMethod = 'childorder'; 
    colormap(ax2, turbo);
    
    for i = 1:nSegments
        for k = 1:nLayers 
            qx = [X_int(i, k), X_int(i+1, k), X_int(i+1, k+1), X_int(i, k+1)];
            qy = [Y_int(i, k), Y_int(i+1, k), Y_int(i+1, k+1), Y_int(i, k+1)];
              
            T_cold = T_matrix(k, i);
            T_hot  = T_matrix(k+1, i);
            
            T_vertices = [T_cold, T_cold, T_hot, T_hot];
            
            patch(ax2, qx, qy, T_vertices, 'EdgeColor', 'none', 'FaceColor', 'interp');
        end
    end
    
    for k = 1:(nLayers + 1)
        if k == 1 || k == (nLayers + 1)
            lw = 1.5; 
            lc = 'k'; 
            ls = '-';
        else
            lw = 1.0; 
            lc = 'k'; 
            ls = '--'; 
        end
        
        plot(ax2, X_int(:, k), Y_int(:, k), 'Color', lc, 'LineStyle', ls, 'LineWidth', lw);
    end
    
    plot(ax2, [min(c_nodes(1,:)), max(c_nodes(1,:))], [0, 0], 'k-.', 'LineWidth', 1);
    
    cb = colorbar(ax2);
    cb.Label.String = 'Layer Temperature (K)';
    cb.Label.FontSize = 11;
    
    xlabel(ax2, 'Axial Position, x (m)');
    ylabel(ax2, 'Radial Position, y (m)');
    
    if all(thicknessScales == 1.0)
        title(ax2, '2D Wall Temperature Distribution (True Scale)');
    elseif isscalar(unique(thicknessScales))
        title(ax2, sprintf('2D Wall Temperature Distribution (Thickness Scaled %gx)', thicknessScales(1)));
    else
        title(ax2, '2D Wall Temperature Distribution (Variable Thickness Scales)');
    end
    xlim(ax2, [min(c_nodes(1,:)), max(c_nodes(1,:))]);

    figure('Name', [prefix 'Mach Number Distribution'], 'Color', 'w', 'WindowStyle', 'docked');
    p = plot(x_mids, M_num, '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    p.MarkerFaceColor = p.Color;
    grid on; box on;
    xlabel('Axial Position, x (m)');
    ylabel('Mach Number, M');
    title('Mach Number Distribution');
    xlim([min(x_mids), max(x_mids)]);

    figure('Name', [prefix 'Convective Heat Transfer'], 'Color', 'w', 'WindowStyle', 'docked');
    p = plot(x_mids, h_conv, '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    p.MarkerFaceColor = p.Color;
    grid on; box on;
    xlabel('Axial Position, x (m)');
    ylabel('h (W/m^2K)');
    title('Bartz Convective Heat Transfer Coefficient');
    xlim([min(x_mids), max(x_mids)]);
    
    figure('Name', [prefix 'Wall Heat Flux'], 'Color', 'w', 'WindowStyle', 'docked');
    p = plot(x_mids, q_flux / 1e6, '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    p.MarkerFaceColor = p.Color;
    grid on; box on;
    xlabel('Axial Position, x (m)');
    ylabel('q'''' (MW/m^2)');
    title('Wall Heat Flux Profile');
    xlim([min(x_mids), max(x_mids)]);

    figure('Name', [prefix 'Area Ratio Distribution'], 'Color', 'w', 'WindowStyle', 'docked');
    p = plot(x_mids, AxAt, '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    p.MarkerFaceColor = p.Color;
    grid on; box on;
    xlabel('Axial Position, x (m)');
    ylabel('A / A_t');
    title('Nozzle Area Ratio Distribution');
    xlim([min(x_mids), max(x_mids)]);

    figure('Name', [prefix 'Gas Temperatures'], 'Color', 'w', 'WindowStyle', 'docked');
    hold on; grid on; box on;
    p1 = plot(x_mids, TotalTemp, '-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Total Temperature (T_t)');
    p1.MarkerFaceColor = p1.Color;
    p2 = plot(x_mids, T_cea, '-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Static Temperature (T_{CEA})');
    p2.MarkerFaceColor = p2.Color;
    xlabel('Axial Position, x (m)');
    ylabel('Temperature (K)');
    title('Core Gas Temperature Distribution');
    legend('Location', 'best');
    xlim([min(x_mids), max(x_mids)]);

    figure('Name', [prefix 'Specific Heat Capacity'], 'Color', 'w', 'WindowStyle', 'docked');
    p = plot(x_mids, cp, '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    p.MarkerFaceColor = p.Color;
    grid on; box on;
    xlabel('Axial Position, x (m)');
    ylabel('c_p (J/kg\cdotK)');
    title('Specific Heat Capacity at Constant Pressure');
    xlim([min(x_mids), max(x_mids)]);

    figure('Name', [prefix 'Specific Heat Ratio'], 'Color', 'w', 'WindowStyle', 'docked');
    p = plot(x_mids, gamma, '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    p.MarkerFaceColor = p.Color;
    grid on; box on;
    xlabel('Axial Position, x (m)');
    ylabel('\gamma');
    title('Ratio of Specific Heats (\gamma)');
    xlim([min(x_mids), max(x_mids)]);

    figure('Name', [prefix 'Prandtl Number'], 'Color', 'w', 'WindowStyle', 'docked');
    p = plot(x_mids, prandtl, '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    p.MarkerFaceColor = p.Color;
    grid on; box on;
    xlabel('Axial Position, x (m)');
    ylabel('Pr');
    title('Prandtl Number Distribution');
    xlim([min(x_mids), max(x_mids)]);

    figure('Name', [prefix 'Dynamic Viscosity'], 'Color', 'w', 'WindowStyle', 'docked');
    p = plot(x_mids, mu, '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    p.MarkerFaceColor = p.Color;
    grid on; box on;
    xlabel('Axial Position, x (m)');
    ylabel('\mu (Pa\cdot s)');
    title('Dynamic Viscosity Distribution');
    xlim([min(x_mids), max(x_mids)]);

    figure('Name', [prefix 'Characteristic Velocity'], 'Color', 'w', 'WindowStyle', 'docked');
    p = plot(x_mids, cstar, '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    p.MarkerFaceColor = p.Color;
    grid on; box on;
    xlabel('Axial Position, x (m)');
    ylabel('c^* (m/s)');
    title('Characteristic Velocity');
    xlim([min(x_mids), max(x_mids)]);

end