% --- Console summary ---
fprintf('\n');
disp(['a                    = ', num2str(a),        ' mm/s/bar^n'])
disp(['a0_uncertainty        = ', num2str(Inc_a),])
disp(['n                    = ', num2str(n),])
disp(['n_uncertainty        = ', num2str(Inc_m),])
disp(['regression R2        = ', num2str(R2),])

fprintf('\n');
fprintf('══════════════════════════════════════════\n');
fprintf('         ENGINE PERFORMANCE SUMMARY       \n');
fprintf('══════════════════════════════════════════\n');
fprintf('  Burn time          : %8.2f  s           \n', time_end_combustion);
fprintf('  Total impulse      : %8.1f  Ns          \n', I_tot_sim);
fprintf('  Design impulse     : %8.1f  Ns          \n', grain.design.I_tot);
fprintf('  Average thrust     : %8.2f  N           \n', T_avg);
fprintf('  Peak thrust        : %8.2f  N  @ %.2fs  \n', T_max, t_Tmax);
fprintf('  Peak chamber press : %8.3f  bar         \n', p_max/1e5);
fprintf('  Design chamber p.  : %8.3f  bar         \n', grain.design.p0/1e5);
grain.SIM_Ispecific = I_tot_sim/(9.80665*grain.formulation.rho_av*(grain.L0*pi*(grain.Re0^2 - grain.Ri0^2)));
fprintf('  Real specific impulse: %8.2f  s           \n', grain.SIM_Ispecific);
fprintf('══════════════════════════════════════════\n\n');

% --- Figure 1: Chamber pressure vs time ---
figure('Name','Chamber Pressure','NumberTitle','off','Color','w');
plot(t_cut, p_cut/1e5, 'b-', 'LineWidth', 1.8);
yline(grain.design.p0/1e5, 'r--', 'LineWidth', 1.2, 'Label', 'Design p_0');
xlabel('Time  [s]',          'FontSize', 12);
ylabel('Chamber pressure  [bar]', 'FontSize', 12);
title('Chamber Pressure vs Time',  'FontSize', 14, 'FontWeight', 'bold');
grid on;  grid minor;
legend('Simulated p_0(t)', 'Design p_0', 'Location', 'best');
set(gca, 'FontSize', 11, 'Box', 'on');

% --- Figure 2: Thrust vs time ---
figure('Name','Thrust Curve','NumberTitle','off','Color','w');
area(t_cut, Thrust_cut/1e3, 'FaceColor', [0.2 0.5 0.8], 'FaceAlpha', 0.25, ...
     'EdgeColor', [0.2 0.5 0.8], 'LineWidth', 1.8);
yline(T_avg/1e3, 'r--', 'LineWidth', 1.2, 'Label', sprintf('T_{avg} = %.1f kN', T_avg/1e3));
yline(grain.design.T_ideal/1e3, 'k:', 'LineWidth', 1.2, ...
      'Label', sprintf('T_{design} = %.1f kN', grain.design.T_ideal/1e3));
xlabel('Time  [s]',    'FontSize', 12);
ylabel('Thrust  [kN]', 'FontSize', 12);
title('Thrust Curve',  'FontSize', 14, 'FontWeight', 'bold');
grid on;  grid minor;
legend('Thrust', 'Average thrust', 'Design thrust', 'Location', 'best');
set(gca, 'FontSize', 11, 'Box', 'on');
