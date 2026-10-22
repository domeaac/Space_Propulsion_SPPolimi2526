function plotWall(Mesh, options)
% PLOTWALL Visualizes the nozzle boundary curve and extruded mesh layers.
%
% Inputs:
%   Mesh                    - Data structure containing layerNodes,
%                             nLayers, and nCurveNodes.
%   options.thicknessScales - (Optional) Scalar or vector specifying the 
%                             magnification factor for layer thicknesses. 
%                             If a vector, length must equal nLayers.
%   options.sectionTitle    - (Optional) String or char scalar to prepend 
%                             to the title of the docked figure tab.

    arguments
        Mesh (1,1) struct
        options.thicknessScales double = []
        options.sectionTitle {mustBeTextScalar} = ""
    end
    
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

    figure('Name', [prefix 'Nozzle Mesh Visualization'], 'Color', 'w', 'WindowStyle', 'docked');
    hold on;
    axis equal;
    grid on;
    xlabel('Axial Coordinate, $x$', 'Interpreter', 'latex');
    ylabel('Radial Coordinate, $y$', 'Interpreter', 'latex');
    
    if all(thicknessScales == 1.0)
        titleStr = 'Structured Nozzle Wall Mesh (True Scale)';
    elseif isscalar(unique(thicknessScales))
        titleStr = sprintf('Structured Nozzle Wall Mesh (Thickness Scaled %gx)', thicknessScales(1));
    else
        titleStr = 'Structured Nozzle Wall Mesh (Variable Thickness Scales)';
    end
    title(titleStr, 'Interpreter', 'latex');

    X_int = zeros(nNodes, nLayers + 1);
    Y_int = zeros(nNodes, nLayers + 1);
    
    l_nodes = Mesh.layerNodes;
    
    % Interface 1 (Outer Base Contour) is extracted directly from page 1
    X_int(:, 1) = l_nodes(1, :, 1)';
    Y_int(:, 1) = l_nodes(2, :, 1)';

    % Compute scaled coordinates by iterating through the self-contained layerNodes pages
    for k = 1:nLayers
        unscaledPrevX = l_nodes(1, :, k)';
        unscaledPrevY = l_nodes(2, :, k)';
        
        unscaledCurrX = l_nodes(1, :, k+1)';
        unscaledCurrY = l_nodes(2, :, k+1)';
        
        dx_layer = unscaledCurrX - unscaledPrevX;
        dy_layer = unscaledCurrY - unscaledPrevY;
        
        X_int(:, k+1) = X_int(:, k) + thicknessScales(k) * dx_layer;
        Y_int(:, k+1) = Y_int(:, k) + thicknessScales(k) * dy_layer;
    end

    % Plot base contour (Cold Outer Wall)
    plot(X_int(:, 1), Y_int(:, 1), 'k.-', 'LineWidth', 1.5, ...
         'MarkerSize', 10, 'DisplayName', 'Outer wall Contour');

    % Plot extruded layer interfaces
    for k = 2:(nLayers + 1)
        plot(X_int(:, k), Y_int(:, k), 'b.-', ...
             'LineWidth', 1, 'MarkerSize', 8, 'HandleVisibility', 'off');
        if k == nLayers+1 % plot hot layer 
             plot(X_int(:, k), Y_int(:, k), 'r.-', 'LineWidth', 1.5, ...
         'MarkerSize', 10, 'DisplayName', 'Inner wall Contour');
        end 


    end

    % Plot transverse structural edges
    for j = 1:nNodes
        plot(X_int(j, :), Y_int(j, :), '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5, ...
             'HandleVisibility', 'off');
    end

    legend('Location', 'best', 'Interpreter', 'latex');
    axis equal; 
    hold off;
end