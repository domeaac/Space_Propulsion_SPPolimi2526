%% Simple MonteCarlo Campaign (Lower-Accuracy Ballistic Model)

if useMontecarlo

% --- Nominal simple-model run (reference for convergence plots) ---
c_star_nom = grain.initial_CEA_run.output.froz.cstar(1) * 1000;   % [mm/s]

dz_nom = @(Ab) ((1e-5)*(a*1e-3)*(grain.formulation.rho_av)*(c_star_nom*1e-3)*Ab./(grain.at*1e6)).^(n/(1-n));
Ab_nom = @(z)  (a*z < (grain.Re0-grain.Ri0)*1e3) .* ...
               (2*pi*(grain.Ri0*1e3 + a*z).*(grain.L0*1e3 - 2*a*z) + ...
                2*pi*((grain.Re0*1e3)^2 - (grain.Ri0*1e3 + a*z).^2));

t0 = 0; tf = 40; dt = 1e-3; z0 = 0; %#ok<NASGU>
f  = @(t, z) dz_nom(Ab_nom(z));
MYRK4_SP_ScriptVersion

dz_num_nom  = dz_nom(Ab_nom(z_vec));
p_sol_nom   = real(dz_num_nom .^ (1/n));
T_sol_nom   = p_sol_nom * 1e5 * grain.at * grain.initial_CEA_run.output.froz.cf(3) * grain.LossesCoeff;
idx_bo_nom  = find(T_sol_nom > 1, 1, 'last');
real_tb_nom = t_vec(idx_bo_nom);
Itot_nom    = trapz(t_vec(1:idx_bo_nom), T_sol_nom(1:idx_bo_nom));   % [Ns]
Tavg_nom    = Itot_nom / real_tb_nom;                                  % [N]
pmax_nom    = max(p_sol_nom(1:idx_bo_nom));                            % [bar]
tburn_nom   = real_tb_nom;                                             % [s]

% --- Monte Carlo setup ---
MC.Tavg  = zeros(1, N_MC^2);
MC.Itot  = zeros(1, N_MC^2);
MC.pmax  = zeros(1, N_MC^2);
MC.tburn = zeros(1, N_MC^2);

% Stochastic properties: only a, n
a_mean  = a;
n_mean  = n;
Sigma_a = Inc_a;
Sigma_n = Inc_n;

% Seeds planting (aleatory uncertainties)
% 1) Create population before starting
MC.a_seeds = a_mean + Sigma_a * randn(1, N_MC);
MC.n_seeds = n_mean + Sigma_n * randn(1, N_MC);

% 2) Create ordered couples
[a_grid, n_grid] = ndgrid(MC.a_seeds, MC.n_seeds);
Couples = [a_grid(:), n_grid(:)];

% 3) Shuffle so seeds are tested in an unordered way to monitor convergence
idx_shuf = randperm(N_MC^2);
Couples_shuffled = Couples(idx_shuf, :);

% Allocate progressive means/deviations for convergence monitoring
MC.Tavgmean  = zeros(1, N_MC^2);
MC.Itotmean  = zeros(1, N_MC^2);
MC.pmaxmean  = zeros(1, N_MC^2);
MC.tburnmean = zeros(1, N_MC^2);
MC.Tavgdev   = zeros(1, N_MC^2);
MC.Itotdev   = zeros(1, N_MC^2);
MC.pmaxdev   = zeros(1, N_MC^2);
MC.tburndev  = zeros(1, N_MC^2);

% USING NOMINAL GEOMETRY COMPUTED WITH NOMINAL a, n (same rationale as
% Montecarlo.m): grain is manufactured at nominal a,n, then burned with
% slightly different ballistic coefficients.

% RK4 time settings (fixed across all iterations)
t0 = 0; tf = 40; dt = 1e-3; z0 = 0;

for k = 1 : N_MC^2

    fprintf("Simple Monte Carlo: iteration %d/%d \n", k, N_MC^2)

    a = Couples_shuffled(k, 1);
    n = Couples_shuffled(k, 2);

    % Simple quasi-steady ballistic model
    dz = @(Ab) ((1e-5)*(a*1e-3)*(grain.formulation.rho_av)*(c_star_nom*1e-3)*Ab./(grain.at*1e6)).^(n/(1-n));
    Ab = @(z)  (a*z < (grain.Re0-grain.Ri0)*1e3) .* ...
               (2*pi*(grain.Ri0*1e3 + a*z).*(grain.L0*1e3 - 2*a*z) + ...
                2*pi*((grain.Re0*1e3)^2 - (grain.Ri0*1e3 + a*z).^2));

    f = @(t, z) dz(Ab(z));

    % Integrate regression depth (MYRK4_SP_ScriptVersion uses t0,tf,dt,f,z0)
    MYRK4_SP_ScriptVersion   % produces t_vec, z_vec; overwrites N internally

    % Post-process: recover p(t) and T(t) from z(t)
    dz_num = dz(Ab(z_vec));
    p_sol  = real(dz_num .^ (1/n));                                              % [bar]
    T_sol  = p_sol * 1e5 * grain.at * grain.initial_CEA_run.output.froz.cf(3) * grain.LossesCoeff; % [N]

    idx_burnout = find(T_sol > 1, 1, 'last');
    if isempty(idx_burnout), idx_burnout = 1; end
    real_tb    = t_vec(idx_burnout);
    real_I_tot = trapz(t_vec(1:idx_burnout), T_sol(1:idx_burnout));   % [Ns]
    T_avg_k    = real_I_tot / real_tb;                                  % [N]
    p_max_k    = max(p_sol(1:idx_burnout));                             % [bar]

    MC.Tavg(k)  = T_avg_k;
    MC.Itot(k)  = real_I_tot;
    MC.pmax(k)  = p_max_k;
    MC.tburn(k) = real_tb;

    % Cumulative means and deviations (Law of Large Numbers convergence)
    MC.Tavgmean(k)  = mean(MC.Tavg(1:k));
    MC.Itotmean(k)  = mean(MC.Itot(1:k));
    MC.pmaxmean(k)  = mean(MC.pmax(1:k));
    MC.tburnmean(k) = mean(MC.tburn(1:k));
    MC.Tavgdev(k)   = std(MC.Tavg(1:k));
    MC.Itotdev(k)   = std(MC.Itot(1:k));
    MC.pmaxdev(k)   = std(MC.pmax(1:k));
    MC.tburndev(k)  = std(MC.tburn(1:k));

end

%% Convergence plots
    figure('Name', 'Simple MC Convergence', 'Color', 'w');

    % --- TOP LEFT: Average Thrust Mean ---
    subplot(4,2,1);
    plot(1:N_MC^2, MC.Tavgmean, 'b', 'LineWidth', 1.5); hold on;
    yline(Tavg_nom, '--r', 'Nominal', 'LineWidth', 1.2);
    grid on; title('Convergence: Mean T_{avg}  (simple model)');
    ylabel('[N]'); xlabel('Iteration');

    % --- TOP RIGHT: Average Thrust Std Dev ---
    subplot(4,2,2);
    plot(1:N_MC^2, MC.Tavgdev, 'k', 'LineWidth', 1.5);
    grid on; title('Convergence: Std Dev T_{avg}  (simple model)');
    ylabel('[N]'); xlabel('Iteration');

    % --- ROW 2 LEFT: Total Impulse Mean ---
    subplot(4,2,3);
    plot(1:N_MC^2, MC.Itotmean, 'b', 'LineWidth', 1.5); hold on;
    yline(Itot_nom, '--r', 'Nominal', 'LineWidth', 1.2);
    grid on; title('Convergence: Mean I_{tot}  (simple model)');
    ylabel('[Ns]'); xlabel('Iteration');

    % --- ROW 2 RIGHT: Total Impulse Std Dev ---
    subplot(4,2,4);
    plot(1:N_MC^2, MC.Itotdev, 'k', 'LineWidth', 1.5);
    grid on; title('Convergence: Std Dev I_{tot}  (simple model)');
    ylabel('[Ns]'); xlabel('Iteration');

    % --- ROW 3 LEFT: Peak Pressure Mean ---
    subplot(4,2,5);
    plot(1:N_MC^2, MC.pmaxmean, 'b', 'LineWidth', 1.5); hold on;
    yline(pmax_nom, '--r', 'Nominal', 'LineWidth', 1.2);
    grid on; title('Convergence: Mean p_{max}  (simple model)');
    ylabel('[bar]'); xlabel('Iteration');

    % --- ROW 3 RIGHT: Peak Pressure Std Dev ---
    subplot(4,2,6);
    plot(1:N_MC^2, MC.pmaxdev, 'k', 'LineWidth', 1.5);
    grid on; title('Convergence: Std Dev p_{max}  (simple model)');
    ylabel('[bar]'); xlabel('Iteration');

    % --- BOTTOM LEFT: Burn Time Mean ---
    subplot(4,2,7);
    plot(1:N_MC^2, MC.tburnmean, 'b', 'LineWidth', 1.5); hold on;
    yline(tburn_nom, '--r', 'Nominal', 'LineWidth', 1.2);
    grid on; title('Convergence: Mean t_{burn}  (simple model)');
    ylabel('[s]'); xlabel('Iteration');

    % --- BOTTOM RIGHT: Burn Time Std Dev ---
    subplot(4,2,8);
    plot(1:N_MC^2, MC.tburndev, 'k', 'LineWidth', 1.5);
    grid on; title('Convergence: Std Dev t_{burn}  (simple model)');
    ylabel('[s]'); xlabel('Iteration');

%% Output Distribution Plots
    Nfull  = N_MC^2;
    mu_T   = MC.Tavgmean(Nfull);  sig_T  = MC.Tavgdev(Nfull);
    mu_I   = MC.Itotmean(Nfull);  sig_I  = MC.Itotdev(Nfull);
    mu_p   = MC.pmaxmean(Nfull);  sig_p  = MC.pmaxdev(Nfull);
    mu_tb  = MC.tburnmean(Nfull); sig_tb = MC.tburndev(Nfull);

    figure('Name','Simple MC – Output Distributions','Color','w');

%    % T_avg
%    subplot(2,2,1);
%    histogram(MC.Tavg, 'Normalization','pdf', ...
%        'FaceColor',[0.20 0.50 0.80],'FaceAlpha',0.5,'EdgeColor','none'); hold on;
%    xv = linspace(mu_T-4*sig_T, mu_T+4*sig_T, 300);
%    plot(xv, normpdf(xv,mu_T,sig_T), 'b-', 'LineWidth', 2);
%    xline(Tavg_nom,        '--r', 'Nominal',   'LineWidth', 1.2);
%    xline(mu_T,            '-b',  '\mu',        'LineWidth', 1.2);
%    xline(mu_T+3*sig_T,   ':k',  '+3\sigma',   'LineWidth', 1.2);
%    xline(mu_T-3*sig_T,   ':k',  '-3\sigma',   'LineWidth', 1.2);
%    xlabel('T_{avg}  [N]'); ylabel('pdf');
%    title('T_{avg} distribution  (simple MC)'); grid on;

    % I_tot
%    subplot(2,2,2);
%    histogram(MC.Itot, 'Normalization','pdf', ...
%        'FaceColor',[0.20 0.70 0.30],'FaceAlpha',0.5,'EdgeColor','none'); hold on;
%    xv = linspace(mu_I-4*sig_I, mu_I+4*sig_I, 300);
%    plot(xv, normpdf(xv,mu_I,sig_I), 'Color',[0 0.6 0], 'LineWidth', 2);
%    xline(Itot_nom,        '--r', 'Nominal',   'LineWidth', 1.2);
%    xline(mu_I,            'Color',[0 0.6 0],  'Label','\mu', 'LineWidth', 1.2);
%    xline(mu_I+3*sig_I,   ':k',  '+3\sigma',   'LineWidth', 1.2);
%    xline(mu_I-3*sig_I,   ':k',  '-3\sigma',   'LineWidth', 1.2);
%    xlabel('I_{tot}  [Ns]'); ylabel('pdf');
%    title('I_{tot} distribution  (simple MC)'); grid on;

    % p_max
    subplot(1,2,1);
    histogram(MC.pmax, 'Normalization','pdf', ...
        'FaceColor',[0.80 0.20 0.20],'FaceAlpha',0.5,'EdgeColor','none'); hold on;
    xv = linspace(mu_p-4*sig_p, mu_p+4*sig_p, 300);
    plot(xv, normpdf(xv,mu_p,sig_p), 'r-', 'LineWidth', 2);
    xline(grain.design.p0/1e5, '--m', 'p_{design}', 'LineWidth', 1.5);
    xline(pmax_nom,            '--r', 'Nominal',     'LineWidth', 1.2);
    xline(mu_p,                '-r',  '\mu',          'LineWidth', 1.2);
    xline(mu_p+3*sig_p,       ':k',  '+3\sigma',     'LineWidth', 1.2);
    xline(mu_p-3*sig_p,       ':k',  '-3\sigma',     'LineWidth', 1.2);
    xlabel('p_{max}  [bar]'); ylabel('pdf');
    title('p_{max} distribution  (simple MC)'); grid on;

    % t_burn
    subplot(1,2,2);
    histogram(MC.tburn, 'Normalization','pdf', ...
        'FaceColor',[0.70 0.40 0.10],'FaceAlpha',0.5,'EdgeColor','none'); hold on;
    xv = linspace(mu_tb-4*sig_tb, mu_tb+4*sig_tb, 300);
    plot(xv, normpdf(xv,mu_tb,sig_tb), 'Color',[0.6 0.3 0], 'LineWidth', 2);
    xline(tburn_nom,        '--r', 'Nominal',  'LineWidth', 1.2);
    xline(mu_tb,            'Color',[0.6 0.3 0], 'Label','\mu', 'LineWidth', 1.2);
    xline(mu_tb+3*sig_tb,  ':k',  '+3\sigma',  'LineWidth', 1.2);
    xline(mu_tb-3*sig_tb,  ':k',  '-3\sigma',  'LineWidth', 1.2);
    xlabel('t_{burn}  [s]'); ylabel('pdf');
    title('t_{burn} distribution  (simple MC)'); grid on;

    sgtitle('Simple MC – Output Distributions', 'FontWeight','bold', 'FontSize', 13);

%% CDF of p_max vs design pressure
    pmax_sorted = sort(MC.pmax);
    cdf_vals    = (1:Nfull) / Nfull;
    P_exceed    = mean(MC.pmax > grain.design.p0/1e5);

    figure('Name','Simple MC – p_{max} CDF','Color','w');
    plot(pmax_sorted, cdf_vals*100, 'r-', 'LineWidth', 2); hold on;
    xline(grain.design.p0/1e5, '--k', 'p_{design}', 'LineWidth', 1.5);
    xline(pmax_nom,            '--b', 'Nominal',     'LineWidth', 1.2);
    yline(50, ':',  '50th pct', 'LineWidth', 1);
    yline(99.73, ':', '99.73% (+3\sigma)', 'LineWidth', 1);
    xlabel('p_{max}  [bar]', 'FontSize', 12);
    ylabel('Cumulative probability  [%]', 'FontSize', 12);
    title(sprintf('CDF of p_{max}  –  P(p_{max} > p_{design}) = %.2f%%', P_exceed*100), ...
          'FontSize', 13, 'FontWeight', 'bold');
    grid on; grid minor;
    set(gca, 'FontSize', 11, 'Box', 'on', 'YLim', [0 100]);

%% Input sensitivity: (a, n) space coloured by p_max
    figure('Name','Simple MC – Input Sensitivity','Color','w');
    scatter(Couples_shuffled(:,1), Couples_shuffled(:,2), 50, MC.pmax, 'filled');
    cb = colorbar;
    cb.Label.String = 'p_{max}  [bar]';
    colormap('jet');
    xline(a_mean, '--w', 'a_{nom}', 'LineWidth', 1.2);
    yline(n_mean, '--w', 'n_{nom}', 'LineWidth', 1.2);
    xlabel('a  [mm s^{-1} bar^{-n}]', 'FontSize', 12);
    ylabel('n  [-]',                  'FontSize', 12);
    title('Input space coloured by p_{max}  (simple MC)', 'FontSize', 13, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 11, 'Box', 'on');

%% 3-sigma console summary
    fprintf('\n');
    fprintf('══════════════════════════════════════════════════════════════\n');
    fprintf('          SIMPLE MC – 3-SIGMA PERFORMANCE ENVELOPE           \n');
    fprintf('══════════════════════════════════════════════════════════════\n');
    fprintf('  T_avg   :  %.1f  ±  %.1f N   →  [%.1f, %.1f] N\n',   mu_T,  3*sig_T,  mu_T-3*sig_T,  mu_T+3*sig_T);
    fprintf('  I_tot   :  %.0f  ±  %.0f Ns  →  [%.0f, %.0f] Ns\n', mu_I,  3*sig_I,  mu_I-3*sig_I,  mu_I+3*sig_I);
    fprintf('  p_max   :  %.3f  ±  %.3f bar →  [%.3f, %.3f] bar\n', mu_p,  3*sig_p,  mu_p-3*sig_p,  mu_p+3*sig_p);
    fprintf('  t_burn  :  %.2f  ±  %.2f s   →  [%.2f, %.2f] s\n',   mu_tb, 3*sig_tb, mu_tb-3*sig_tb, mu_tb+3*sig_tb);
    fprintf('  P(p_max > p_design) = %.2f%%\n', P_exceed*100);
    fprintf('══════════════════════════════════════════════════════════════\n\n');

end
