%% MonteCarlo Campaign

if useMontecarlo 

Tavg_nom  = T_avg;
Itot_nom  = I_tot_sim;
pmax_nom  = p_max;
tburn_nom = time_end_combustion;

N = 30;
MC.Tavg  = zeros(1,N^2);
MC.Itot  = zeros(1,N^2);
MC.pmax  = zeros(1,N^2);
MC.tburn = zeros(1,N^2);
% MC.tburn = zeros(1,N); already implied by Tavg,Itot

% Stochastic properties
% only a, n
a_mean = a; % Sample mean estimation (Stima della media campionaria)
n_mean = n; % Sample mean estimation
Sigma_a = Inc_a; % Sample deviation estimation
Sigma_n = Inc_n; % Sample deviation estimation

% Seeds planting
%  aleatory uncertaintes (production, non-epistemic)
% because natural variability of the problem.

% Order of plating is importance to obtain method convergence

% 1) Create population before starting
MC.a_seeds =  a_mean + Sigma_a * randn(1, N);
MC.n_seeds =  n_mean + Sigma_n * randn(1, N);

% 2) Create ordered couples of seeds
[a_grid, n_grid] = ndgrid(MC.a_seeds, MC.n_seeds);
Couples = [a_grid(:), n_grid(:)]; % 2xN

% 3) Feed seeds into shuffling functions: shuffle the data son seeds are tested in an unordered way
% to properly monitor Monte Carlo convergence
idx = randperm(N^2);
Couples_shuffled = Couples(idx,:);


% Allocate progressive means for Monte Carlo convergence
    MC.Tavgmean  = zeros(1,N^2);
    MC.Itotmean  = zeros(1,N^2);
    MC.pmaxmean  = zeros(1,N^2);
    MC.tburnmean = zeros(1,N^2);

    MC.Tavgdev  = zeros(1,N^2);
    MC.Itotdev  = zeros(1,N^2);
    MC.pmaxdev  = zeros(1,N^2);
    MC.tburndev = zeros(1,N^2);

  
% USING NOMINAL GEOMETRY COMPUTED WITH NOMINAL a, n, NOT RECOMPUTING
% GEOMETRY FOR EVERY a, n: indeed realistically you then build it with that
% geometry computed using nominal a,n, but then a, n are slightly
% different, not the geometry (if perfect manufacturing).
% Scope: You have computed the nominal engine and you have the nominal
% geometry; you manufacture the geometry which is corresponding to the
% nominal a,n, and then you burn it, but when you burn it a n are slightly
% different, while geometry is fixed, a,n now do not influence it like in
% the previous lines. And you want to see, in this real bench case
% scenario, what average values of T_av, Itot, delta T you will get (and
% also dev), which are not nec the ones nominal! Because a,n, do not
% influence the grain. Those values are the ones the Montecarlo progressive
% means converge to, not nec the nominal values. Ok.

    for k = 1 : N^2 % Run ballistic code

        fprintf("Monte Carlo: iteration %d/%d \n",k,N^2)

        a = Couples_shuffled(k,1);
        n = Couples_shuffled(k,2);

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
    plot(1:N^2, MC.Tavgmean, 'b', 'LineWidth', 1.5); hold on;
    yline(Tavg_nom, '--r', 'Nominal', 'LineWidth', 1.2);
    grid on; title('Convergence: Mean T_{avg}');
    ylabel('[N]'); xlabel('Iteration');

    % --- TOP RIGHT: Average Thrust Std Dev ---
    subplot(4,2,2);
    plot(1:N^2, MC.Tavgdev, 'k', 'LineWidth', 1.5);
    grid on; title('Convergence: Std Dev T_{avg}');
    ylabel('[N]'); xlabel('Iteration');

    % --- ROW 2 LEFT: Total Impulse Mean ---
    subplot(4,2,3);
    plot(1:N^2, MC.Itotmean, 'b', 'LineWidth', 1.5); hold on;
    yline(Itot_nom, '--r', 'Nominal', 'LineWidth', 1.2);
    grid on; title('Convergence: Mean I_{tot}');
    ylabel('[Ns]'); xlabel('Iteration');

    % --- ROW 2 RIGHT: Total Impulse Std Dev ---
    subplot(4,2,4);
    plot(1:N^2, MC.Itotdev, 'k', 'LineWidth', 1.5);
    grid on; title('Convergence: Std Dev I_{tot}');
    ylabel('[Ns]'); xlabel('Iteration');

    % --- ROW 3 LEFT: Peak Pressure Mean ---
    subplot(4,2,5);
    plot(1:N^2, MC.pmaxmean, 'b', 'LineWidth', 1.5); hold on;
    yline(pmax_nom, '--r', 'Nominal', 'LineWidth', 1.2);
    grid on; title('Convergence: Mean p_{max}');
    ylabel('[Pa]'); xlabel('Iteration');

    % --- ROW 3 RIGHT: Peak Pressure Std Dev ---
    subplot(4,2,6);
    plot(1:N^2, MC.pmaxdev, 'k', 'LineWidth', 1.5);
    grid on; title('Convergence: Std Dev p_{max}');
    ylabel('[Pa]'); xlabel('Iteration');

    % --- BOTTOM LEFT: Burn Time Mean ---
    subplot(4,2,7);
    plot(1:N^2, MC.tburnmean, 'b', 'LineWidth', 1.5); hold on;
    yline(tburn_nom, '--r', 'Nominal', 'LineWidth', 1.2);
    grid on; title('Convergence: Mean t_{burn}');
    ylabel('[s]'); xlabel('Iteration');

    % --- BOTTOM RIGHT: Burn Time Std Dev ---
    subplot(4,2,8);
    plot(1:N^2, MC.tburndev, 'k', 'LineWidth', 1.5);
    grid on; title('Convergence: Std Dev t_{burn}');
    ylabel('[s]'); xlabel('Iteration');

end
% In media quindi il motore performerà diversamente dal nominale, perchè
% sto imponendo proprietà ballistiche variabili ma che non influiscono
% sulla geometria, diversamente da prima, ora H, etc. sono fissati, 
% per una perfomrance media minore/maggiore MC.Tavgmean = , etc.
% quindi quando lo brucerai otterrai questi valori con media e
% deviazionestd... quindi hai la prb del 99.97% di trovare una T =
% Tavg+3sigma, t = etc. considerando pero solo incertezza su a, n in realta
% anche altre fonti e in realtà modello non esatto quindi gli esiti dei
% burn daranno valori con media e deviazione leggermente diversi sicuro.
% Inoltre avere piu dati su a,n permette una mu,sigma di a n minore e
% quindi migliore predizione degli esiti dei burn (scrivi negli improvements).
% SALVA VALORI
% SALVA PLOTS