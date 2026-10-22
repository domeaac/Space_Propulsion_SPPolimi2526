function plot_engine_xsection(L0, Ri0, Re0, At, Ae, Lconv, Ldiv, nozzle_conv_angle_degs, nozzle_div_angle_degs)
    % plot_engine_xsection - Cross-section diagram of full engine
    %
    % Inputs:
    %   L0, Ri0, Re0    - Grain dimensions [mm]
    %   At, Ae          - Throat and exit areas [mm^2]
    %   Lconv, Ldiv     - Convergent and divergent nozzle lengths [mm]
    %   nozzle_conv_angle_degs, nozzle_div_angle_degs - Nozzle angles [deg]
    
    Dt = 2*sqrt(At/pi);
    De = 2*sqrt(Ae/pi);
    
    figure('Position', [100, 100, 1400, 600]);
    hold on; axis equal; grid on;
    
    % ===== GRAIN OUTER CASING =====
    x_casing = [0, L0, L0, 0];
    y_casing = [Re0, Re0, -Re0, -Re0];
    patch(x_casing, y_casing, [0.8 0.8 0.8], 'EdgeColor', 'k', 'LineWidth', 2, 'FaceAlpha', 0.3);
    
    % ===== PROPELLANT (outer grain) =====
    x_prop = [0, L0, L0, 0];
    y_prop = [Re0, Re0, -Re0, -Re0];
    patch(x_prop, y_prop, [1 0.7 0.4], 'EdgeColor', [0.5 0.3 0], 'LineWidth', 2, 'FaceAlpha', 0.8);
    
    % ===== CORE (burned cavity) =====
    x_core = [0, L0, L0, 0];
    y_core = [Ri0, Ri0, -Ri0, -Ri0];
    patch(x_core, y_core, 'white', 'EdgeColor', 'r', 'LineWidth', 1.5);
    
    % ===== CONVERGENT NOZZLE =====
    x_conv = [L0, L0 + Lconv];
    y_conv_top = linspace(Re0, Dt/2, 50);
    x_conv_top = linspace(L0, L0 + Lconv, 50);
    y_conv_bot = -y_conv_top;
    x_conv_bot = x_conv_top;
    
    patch([x_conv_top, fliplr(x_conv_bot)], [y_conv_top, fliplr(y_conv_bot)], ...
        [0.7 0.7 0.9], 'EdgeColor', 'k', 'LineWidth', 2, 'FaceAlpha', 0.6);
    
    % ===== DIVERGENT NOZZLE =====
    x_div = [L0 + Lconv, L0 + Lconv + Ldiv];
    y_div_top = linspace(Dt/2, De/2, 50);
    x_div_top = linspace(L0 + Lconv, L0 + Lconv + Ldiv, 50);
    y_div_bot = -y_div_top;
    x_div_bot = x_div_top;
    
    patch([x_div_top, fliplr(x_div_bot)], [y_div_top, fliplr(y_div_bot)], ...
        [0.5 0.5 0.7], 'EdgeColor', 'k', 'LineWidth', 2, 'FaceAlpha', 0.6);
    
    % ===== ANNOTATIONS =====
    % Grain
    text(L0/2, Re0 + 30, sprintf('GRAIN: L=%.0f mm', L0), ...
        'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
    text(L0/2, -Ri0/2, sprintf('Core R_i=%.1f mm', Ri0), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'red', 'FontWeight', 'bold');
    
    % Nozzle
    text(L0 + Lconv/2, Dt/2 + 40, ...
        sprintf('CONV\n\\theta=%.0f°\nL=%.1f mm', nozzle_conv_angle_degs, Lconv), ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
    text(L0 + Lconv + Ldiv/2, De/2 + 50, ...
        sprintf('DIV\n\\theta=%.0f°\nL=%.1f mm\nA_e/A_t=%.2f', ...
        nozzle_div_angle_degs, Ldiv, Ae/At), ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
    
    % Dimensions
    text(L0 + Lconv + Ldiv + 100, 0, sprintf('D_e=%.1f mm', De), ...
        'FontSize', 10, 'FontWeight', 'bold');
    text(L0 + Lconv/2, -Dt/2 - 40, sprintf('D_t=%.1f mm', Dt), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
    
    xlabel('Axial Distance [mm]', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Radial Distance [mm]', 'FontSize', 12, 'FontWeight', 'bold');
    title('Engine Cross-Section (Axisymmetric)', 'FontSize', 14, 'FontWeight', 'bold');
    
    xlim([-50, L0 + Lconv + Ldiv + 150]);
    ylim([-(Re0 + 100), Re0 + 100]);
    set(gca, 'FontSize', 11);
end