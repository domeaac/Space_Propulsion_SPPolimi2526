%% MOTOR Analysis
plots = true;

c_star = grain.initial_CEA_run.output.froz.cstar(1)*1000;   %mm/s

% CONSERVATION OF MASS ('Pressure Equation') 
dz = @(Ab) ((1e-5)*(a*1e-3)*(grain.formulation.rho_av)*(c_star*1e-3)*Ab./(grain.at*1e6)).^(n/(1-n)); % dz = p^n
% GRAIN GEOMETRY 
Ab = @(z) (a*z<(grain.Re0-grain.Ri0)*1e3).*(2*pi*(grain.Ri0*1e3 + a*z).*(grain.L0*1e3 - 2*a*z) + 2*pi*((grain.Re0*1e3)^2 - (grain.Ri0*1e3 + a*z).^2));
Ab0 = Ab(0);

% I.C. for Ab(t):
% Ab0 = 2*pi*Ri0*L0 + 2*pi*(Re0^2 - Ri0^2);

% time data
t0 = 0;          % [s]
tf = 40;         % [s]
dt = 1e-3;       % passo temporale [s]
f = @(t, z) dz(Ab(z));   
z0 = 0;          % condizione iniziale [s]

MYRK4_SP_ScriptVersion

% Post-processing: ricava p(t) da z(t)
dz_num = dz(Ab(z_vec));
p_sol  = real(dz_num .^ (1/n));

T_sol = p_sol*1e5*grain.at*grain.initial_CEA_run.output.froz.cf(3)*grain.LossesCoeff;
real_I_tot = trapz(T_sol)*dt;

[p_max, idx_max] = max(p_sol);
[T_max, idx_T_max] = max(T_sol);

% Extract actual burn time
% disp(max(T_sol))   % check if thrust ever exceeds 1
% disp(isreal(T_sol)) % check if real

idx_burnout = find(T_sol > 1, 1, 'last');  % last meaningful thrust
real_tb = t_vec(idx_burnout);

% Compute from integration
real_I_tot = trapz(t_vec(1:idx_burnout), T_sol(1:idx_burnout));
T_avg = real_I_tot / real_tb;

% Propellant mass - unsteady solver
% The solver runs to full burnout (grain.x/y are a local copy inside the function),
% so all initial propellant is consumed.
V_prop_init = pi * (grain.Re0^2 - grain.Ri0^2) * grain.L0;   % [m^3] initial grain volume
m_prop_adv  = grain.formulation.rho_av * V_prop_init;         % [kg]
Isp_adv     = I_tot_sim / (m_prop_adv * 9.81);                % [s]

% Propellant mass - quasi-steady solver
% z_vec integrates dz/dt = p^n, so a*z [mm] is the web burned at each instant.
e_sim_m    = a * z_vec(idx_burnout) * 1e-3;                   % [m] web burned at burnout
Ri_f       = grain.Ri0 + e_sim_m;                             % [m] final inner radius
L_f        = max(grain.L0 - 2*e_sim_m, 0);                    % [m] final grain length
V_prop_rem = pi * (grain.Re0^2 - Ri_f^2) * L_f;              % [m^3] unburned volume
m_prop_sim = grain.formulation.rho_av * (V_prop_init - V_prop_rem);  % [kg]
Isp_sim    = real_I_tot / (m_prop_adv * 9.81);                % [s]

% CEA frozen-flow reference Isp (vacuum, at nozzle exit plane, with divergence loss)
Isp_CEA    = grain.design.Isp; % [s]

%% SIMULATION COMPARISON SUMMARY
fprintf('\n');
fprintf('══════════════════════════════════════════════════════════════════\n');
fprintf('                 SIMULATION COMPARISON SUMMARY                    \n');
fprintf('══════════════════════════════════════════════════════════════════\n');
fprintf('  CC dynamics model       Unsteady solver        Quasi-steady solver     \n');
fprintf('  ──────────────────────────────────────────────────────────────  \n');
fprintf('  Burn time          :    %8.2f s          %8.2f s        \n', time_end_combustion, real_tb);
fprintf('  Total impulse      :    %8.1f Ns         %8.1f Ns       \n', I_tot_sim, real_I_tot);
fprintf('  Average thrust     :    %8.2f N          %8.2f N        \n', I_tot_sim/time_end_combustion, T_avg);
fprintf('  Peak thrust        :    %8.2f N          %8.2f N        \n', max(Thrust_cut), T_max);
fprintf('  Peak chamber press :    %8.3f bar        %8.3f bar      \n', max(p_cut)/1e5, p_max);
fprintf('  Design chamber p.  :    %8.3f bar        %8.3f bar      \n', grain.design.p0/1e5, grain.design.p0/1e5);
fprintf('  ──────────────────────────────────────────────────────────────  \n');
fprintf('  Specific impulse   :    %8.2f s           %8.2f s       \n', Isp_adv, Isp_sim);
fprintf('  Isp  (CEA ref, no 2D losses):    %8.2f s                          \n', Isp_CEA);
fprintf('  Propellant mass    :    %8.3f kg          %8.3f kg      \n', m_prop_adv, m_prop_sim);
fprintf('══════════════════════════════════════════════════════════════════\n\n');

if plots
    col_adv  = [0.00 0.35 0.75];   % deep blue  – Unsteady solver curves
    col_adv2 = [0.40 0.70 1.00];   % light blue – Unsteady solver annotation lines
    col_sim  = [0.80 0.10 0.10];   % deep red   – Quasi-steady solver curves
    col_sim2 = [1.00 0.45 0.45];   % light red  – Quasi-steady solver annotation lines

    % Pre-compute indices for advanced model peaks
    [~, idx_pmax_adv]  = max(p_cut);
    [~, idx_Tmax_adv]  = max(Thrust_cut);

    % --- Pressure Trace ---
    figure('Name','Pressure Trace Comparison','NumberTitle','off','Color','w');
    hp1 = plot(t_cut, p_cut/1e5, '-',  'Color', col_adv, 'LineWidth', 1.8); hold on;
    hp2 = plot(t_vec(1:idx_burnout), p_sol(1:idx_burnout), '--', 'Color', col_sim, 'LineWidth', 1.8);
    hp3 = yline(grain.design.p0/1e5, 'k:', 'LineWidth', 1.2, 'Label', 'Design p_0');
    % Burnout points
    hp4 = plot(time_end_combustion, p_cut(end)/1e5,        's', 'Color', col_adv2, 'MarkerFaceColor', col_adv2, 'MarkerSize', 8);
    hp5 = plot(real_tb,             p_sol(idx_burnout),    's', 'Color', col_sim2, 'MarkerFaceColor', col_sim2, 'MarkerSize', 8);
    % Pmax points
    hp6 = plot(t_cut(idx_pmax_adv), p_cut(idx_pmax_adv)/1e5, '^', 'Color', col_adv2, 'MarkerFaceColor', col_adv2, 'MarkerSize', 8);
    hp7 = plot(t_vec(idx_max),      p_max,                    '^', 'Color', col_sim2, 'MarkerFaceColor', col_sim2, 'MarkerSize', 8);
    xlabel('Time  [s]', 'FontSize', 12);
    ylabel('Chamber pressure  [bar]', 'FontSize', 12);
    title('Chamber Pressure Comparison', 'FontSize', 14, 'FontWeight', 'bold');
    grid on; grid minor;
    legend([hp1 hp2 hp3 hp4 hp5 hp6 hp7], ...
           'Unsteady solver', 'Quasi-steady solver', 'Design p_0', ...
           'Burnout (adv.)', 'Burnout (simple)', 'p_{max} (adv.)', 'p_{max} (simple)', ...
           'Location', 'best');
    set(gca, 'FontSize', 11, 'Box', 'on');

    % --- Thrust Curve ---
    figure('Name','Thrust Comparison','NumberTitle','off','Color','w');
    ht1 = plot(t_cut, Thrust_cut, '-',  'Color', col_adv, 'LineWidth', 1.8); hold on;
    ht2 = plot(t_vec(1:idx_burnout), T_sol(1:idx_burnout), '--', 'Color', col_sim, 'LineWidth', 1.8);
    % Burnout points
    ht3 = plot(time_end_combustion, Thrust_cut(end),    's', 'Color', col_adv2, 'MarkerFaceColor', col_adv2, 'MarkerSize', 8);
    ht4 = plot(real_tb,             T_sol(idx_burnout), 's', 'Color', col_sim2, 'MarkerFaceColor', col_sim2, 'MarkerSize', 8);
    % Tmax points
    ht5 = plot(t_cut(idx_Tmax_adv), Thrust_cut(idx_Tmax_adv), '^', 'Color', col_adv2, 'MarkerFaceColor', col_adv2, 'MarkerSize', 8);
    ht6 = plot(t_vec(idx_T_max),    T_max,                     '^', 'Color', col_sim2, 'MarkerFaceColor', col_sim2, 'MarkerSize', 8);
    % Taverage points (placed at mid-burn time)
    ht7 = plot(time_end_combustion/2, I_tot_sim/time_end_combustion, 'd', 'Color', col_adv2, 'MarkerFaceColor', col_adv2, 'MarkerSize', 8);
    ht8 = plot(real_tb/2,             T_avg,                         'd', 'Color', col_sim2, 'MarkerFaceColor', col_sim2, 'MarkerSize', 8);
    xlabel('Time  [s]', 'FontSize', 12);
    ylabel('Thrust  [N]', 'FontSize', 12);
    title('Thrust Curve Comparison', 'FontSize', 14, 'FontWeight', 'bold');
    grid on; grid minor;
    legend([ht1 ht2 ht3 ht4 ht5 ht6 ht7 ht8], ...
           'Unsteady solver', 'Quasi-steady solver', ...
           'Burnout (adv.)', 'Burnout (simple)', ...
           'T_{max} (adv.)', 'T_{max} (simple)', ...
           'T_{avg} (adv.)', 'T_{avg} (simple)', ...
           'Location', 'best');
    set(gca, 'FontSize', 11, 'Box', 'on');

end






