clear;
load("montecarlo2004DomeWorkspace.mat")

% Convergence plots
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
