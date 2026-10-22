%% MonteCarlo Campaign

if useMontecarlo 

Tavg_nom  = T_avg;
Itot_nom  = I_tot_sim;
pmax_nom  = p_max;
tburn_nom = time_end_combustion;

N_MC = 900;
MC.Tavg  = zeros(1,N_MC);
MC.Itot  = zeros(1,N_MC);
MC.pmax  = zeros(1,N_MC);
MC.tburn = zeros(1,N_MC);

% Stochastic properties — only a, n
a_mean = a;
n_mean = n;
Sigma_a = Inc_a;
Sigma_n = Inc_n;


% Allocate progressive means for Monte Carlo convergence
    MC.Tavgmean  = zeros(1,N_MC);
    MC.Itotmean  = zeros(1,N_MC);
    MC.pmaxmean  = zeros(1,N_MC);
    MC.tburnmean = zeros(1,N_MC);

    MC.Tavgdev  = zeros(1,N_MC);
    MC.Itotdev  = zeros(1,N_MC);
    MC.pmaxdev  = zeros(1,N_MC);
    MC.tburndev = zeros(1,N_MC);

    for k = 1 : N_MC % Run ballistic code

        fprintf("Monte Carlo: iteration %d/%d \n",k,N_MC)

        a = a_mean + Sigma_a * randn();
        n = n_mean + Sigma_n * randn();

        grain.r_b = @(p) a*p.^n;

   % Run ballistics (do not print processing: use evalc)
   [~, t, sol, mass_rate, Thrust, time_end_combustion] = evalc('reduced_RK_2_frozen_Z(grain, grain.pe + 10e5, timestep, tend)'); 
    % post-processing of results:
    % --- Trim to combustion window ---
        idx_end = round(time_end_combustion / timestep);
        Thrust_cut = Thrust(1:idx_end);
        t_cut      = t(1:idx_end)  ;    
        p_cut      = sol(1:idx_end);  % chamber pressure [Pa]

    % --- Performance metrics ---
        I_tot_sim = trapz(t_cut, Thrust_cut);                 % [Ns]  total impulse
        T_avg     = I_tot_sim / time_end_combustion;          % [N]   average thrust
        [T_max, idx_Tmax] = max(Thrust_cut);                  % [N]   peak thrust
        p_max     = max(p_cut);                               % [Pa]  peak pressure
        t_Tmax    = t_cut(idx_Tmax);                          % [s]

        MC.Tavg (k)  = T_avg;
        MC.Itot (k)  = I_tot_sim;
        MC.pmax (k)  = p_max;
        MC.tburn(k)  = time_end_combustion;

        % Monte Carlo convergence

        % Cumulative means (-> mean for Legge Grandi Num)
        MC.Tavgmean(k)  = mean(MC.Tavg(1:k));
        MC.Itotmean(k)  = mean(MC.Itot(1:k));
        MC.pmaxmean(k)  = mean(MC.pmax(1:k));
        MC.tburnmean(k) = mean(MC.tburn(1:k));

        % Cumulative devations (-> mean for Legge Grandi Num)
        MC.Tavgdev(k)  = std(MC.Tavg(1:k));
        MC.Itotdev(k)  = std(MC.Itot(1:k));
        MC.pmaxdev(k)  = std(MC.pmax(1:k));
        MC.tburndev(k) = std(MC.tburn(1:k));

    end

%%
    figure('Name', 'Monte Carlo Convergence', 'Color', 'w');

    % --- TOP LEFT: Average Thrust Mean ---
    subplot(4,2,1);
    plot(1:N_MC, MC.Tavgmean, 'b', 'LineWidth', 1.5); hold on;
    yline(Tavg_nom, '--r', 'Nominal', 'LineWidth', 1.2);
    grid on; title('Convergence: Mean T_{avg}');
    ylabel('[N]'); xlabel('Iteration');

    % --- TOP RIGHT: Average Thrust Std Dev ---
    subplot(4,2,2);
    plot(1:N_MC, MC.Tavgdev, 'k', 'LineWidth', 1.5);
    grid on; title('Convergence: Std Dev T_{avg}');
    ylabel('[N]'); xlabel('Iteration');

    % --- ROW 2 LEFT: Total Impulse Mean ---
    subplot(4,2,3);
    plot(1:N_MC, MC.Itotmean, 'b', 'LineWidth', 1.5); hold on;
    yline(Itot_nom, '--r', 'Nominal', 'LineWidth', 1.2);
    grid on; title('Convergence: Mean I_{tot}');
    ylabel('[Ns]'); xlabel('Iteration');

    % --- ROW 2 RIGHT: Total Impulse Std Dev ---
    subplot(4,2,4);
    plot(1:N_MC, MC.Itotdev, 'k', 'LineWidth', 1.5);
    grid on; title('Convergence: Std Dev I_{tot}');
    ylabel('[Ns]'); xlabel('Iteration');

    % --- ROW 3 LEFT: Peak Pressure Mean ---
    subplot(4,2,5);
    plot(1:N_MC, MC.pmaxmean, 'b', 'LineWidth', 1.5); hold on;
    yline(pmax_nom, '--r', 'Nominal', 'LineWidth', 1.2);
    grid on; title('Convergence: Mean p_{max}');
    ylabel('[Pa]'); xlabel('Iteration');

    % --- ROW 3 RIGHT: Peak Pressure Std Dev ---
    subplot(4,2,6);
    plot(1:N_MC, MC.pmaxdev, 'k', 'LineWidth', 1.5);
    grid on; title('Convergence: Std Dev p_{max}');
    ylabel('[Pa]'); xlabel('Iteration');

    % --- BOTTOM LEFT: Burn Time Mean ---
    subplot(4,2,7);
    plot(1:N_MC, MC.tburnmean, 'b', 'LineWidth', 1.5); hold on;
    yline(tburn_nom, '--r', 'Nominal', 'LineWidth', 1.2);
    grid on; title('Convergence: Mean t_{burn}');
    ylabel('[s]'); xlabel('Iteration');

    % --- BOTTOM RIGHT: Burn Time Std Dev ---
    subplot(4,2,8);
    plot(1:N_MC, MC.tburndev, 'k', 'LineWidth', 1.5);
    grid on; title('Convergence: Std Dev t_{burn}');
    ylabel('[s]'); xlabel('Iteration');

end
