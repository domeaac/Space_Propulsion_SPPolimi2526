%% TARGETS
I_tot_target = 2.5e6;   % [Ns]
T_avg_target = 100e3;   % [N]


%% FMINCON  — x = [Re0 (m), L0 (m), At (m²), Ri0 (m)]
x0 = [grain.Re0, grain.L0, grain.at, grain.Ri0];
lb = [0.10, 0.20, 1e-4, 0.05];
ub = [1.00, 3.00, 0.10, 0.80];

% nonlinear constraint: Ri0 < Re0
nlcon = @(x) deal(x(4) - x(1), []);

opts = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
    'OptimalityTolerance', 1e-4, 'StepTolerance', 1e-6, 'MaxFunctionEvaluations', 300);

[x_opt, fval, exitflag] = fmincon( ...
    @(x) obj(x, grain, I_tot_target, T_avg_target), ...
    x0, [], [], [], [], lb, ub, nlcon, opts);

fprintf('\nExit: %d  |  J = %.4e  |  Re0=%.3fm  L0=%.3fm  At=%.4fm2  Ri0=%.3fm\n', ...
    exitflag, fval, x_opt(1), x_opt(2), x_opt(3), x_opt(4));


%% OBJECTIVE FUNCTION
function J = obj(x, grain, I_tot_target, T_avg_target)
    grain.Re0 = x(1);
    grain.L0  = x(2);
    grain.at  = x(3);
    grain.Ri0 = x(4);
    grain.ae = grain.at*grain.aeat;

    grain = initGrainMesh(grain);

    try
        [t, ~, ~, Thrust, t_end] = reduced_RK_2_frozen_Z(grain, grain.design.p0, 0.1, 200);
    catch % safeguard the optimizer !!!
        J = 1e6;  return
    end

    if t_end == 0 || isempty(Thrust);  J = 1e6;  return;  end

    idx   = max(1, round(t_end / 0.1));
    I_sim = trapz(t(1:idx), Thrust(1:idx));
    T_avg = I_sim / t_end;

    J = 10*((I_sim - I_tot_target) / I_tot_target)^2 + ((T_avg - T_avg_target) / T_avg_target)^2;
    fprintf('  Re0=%.3f L0=%.3f At=%.4f Ri0=%.3f  =>  I=%7.0fNs  Tavg=%6.1fN  J=%.3e\n', ...
        x(1), x(2), x(3), x(4), I_sim, T_avg, J);
end
