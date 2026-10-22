clear, clc, close all
restoredefaultpath
addpath(genpath(pwd))
plots = true;

%% DESIGN SPACE VATIABLES

% temp
grain.temperature = 298.15;                  % [K]

% pressures
grain.design.p0 = 70*1e5;                           % [Pa]
grain.pe = 101325;                                  % [Pa]
grain.design.p0_over_pe = grain.design.p0/grain.pe; % []
 
% performance
grain.design.I_tot = 2.5*1e6;                      % [Ns] 
grain.design.T_ideal = 100*1e3;                    % [N]

% propellant data
grain.formulation.AP_mass_fraction = 80;   % []
grain.formulation.HTPB_mass_fraction = 20; % []

% geometry
grain.nozzle_div_angle_degs = 15;           
grain.nozzle_conv_angle_degs = 45;


%% CEA Analysis of propellant formulation  - Variable AP and HTPB percentages

% run CEA to init grain design parameters
grain.initial_CEA_run = wrp_rkt_AP_HTPB_imposing_pipe(grain.design.p0/1e5, grain.design.p0_over_pe, grain.temperature, ...
    grain.formulation.AP_mass_fraction, grain.formulation.HTPB_mass_fraction);

% save parameters
grain.aeat = grain.initial_CEA_run.output.froz.aeat(3);
grain.cstar = grain.initial_CEA_run.output.froz.cstar(1);
grain.design.gamma_chamber = grain.initial_CEA_run.output.froz.gamma(1);  

%% PROVIDED DATA ANALYSIS

% Experimental data loading
load("provided_data_file.mat")
p_val = provided_data(:,1);
r_val = provided_data(:,2);

% Compute Vielle Law Coeffiecients
UncertaintyGiorgio
disp(['a                    = ', num2str(a),        ' mm/s/bar^n'])
disp(['n                    = ', num2str(n),])
disp(['n_uncertainty        = ', num2str(Inc_m),])
disp(['regression R2        = ', num2str(R2),])

% VIEILLE LAW 
grain.r_b = @(p) a*p.^n;

%% BATES MOTOR Design

% process design data
grain.LossesCoeff = (1+cos(deg2rad(grain.nozzle_div_angle_degs)))/2; % []
grain.design.T = grain.design.T_ideal/grain.LossesCoeff;                     % [N]
grain.design.tb = grain.design.I_tot / grain.design.T_ideal;                                    % [s] 

% Mass Properties:
grain.formulation.mperc_vec = [grain.formulation.AP_mass_fraction, grain.formulation.HTPB_mass_fraction]/100;  % [] OX, FU percentages
grain.formulation.rho_vec   = [1950, 920];        % [kg/m^3]
grain.formulation.rho_av = 1/sum(grain.formulation.mperc_vec./grain.formulation.rho_vec);  % [kg/m^3] average propellant density

% Burning rate in nominal conditions;
grain.design.rb = grain.r_b(grain.design.p0/1e5);                      % [mm/s]

%% BERA code integration - Bates RM Design
grain.design.ct = grain.initial_CEA_run.output.froz.cf(3);
grain.design.m_dot = grain.design.T /(grain.design.ct*grain.cstar);         % [kg/s]
% compute ideal At
grain.at = grain.cstar * grain.design.m_dot / grain.design.p0; % [m^2]

% initial geometry solver use
initial_conditions = [1.5 , 0.5 , 0.8];    % L0, Ri0, Re0, [m]
Bates = solveBates(grain, grain.design.rb/1e3, grain.design.m_dot, grain.design.tb, grain.formulation.rho_av, "x0", initial_conditions); 

% save results
grain.L0  = Bates.h;   
grain.Ri0 = Bates.D_int/2;
grain.Re0 = Bates.D_ext/2;


%% PERFORMANCE SIMULATION
% RESULTS FROM OPTIMIZATION 
 % grain.L0 = 1.640;
 % grain.Re0 = 0.451;
 % grain.Ri0 = 0.2582;
 % grain.at = 0.0107;

clc
grain = initGrainMesh(grain);
timestep = 0.1;
tend = 100;
% complete model - 4*x runtime
% [t,sol, mass_rate, Thrust, time_end_combustion] = my_RK_2_frozen_Z(grain, grain.design.pe + 10e5, timestep, tend);
% constant CEA output in RK  - x runtime
[t,sol, mass_rate, Thrust, time_end_combustion] = reduced_RK_2_frozen_Z(grain, grain.pe + 10e5, timestep, tend);

%% POSTPROCESS RESULTS

% --- Trim to combustion window ---
idx_end = round(time_end_combustion / timestep);
Thrust_cut = Thrust(1:idx_end);
t_cut      = t(1:idx_end);
p_cut      = sol(1:idx_end);           % chamber pressure [Pa]

% --- Performance metrics ---
I_tot_sim = trapz(t_cut, Thrust_cut);                % [Ns]  total impulse
T_avg     = I_tot_sim / time_end_combustion;          % [N]   average thrust
[T_max, idx_Tmax] = max(Thrust_cut);                  % [N]   peak thrust
p_max     = max(p_cut);                               % [Pa]  peak pressure
t_Tmax    = t_cut(idx_Tmax);                          % [s]

if plots
PlotAndPrint
end
%% NOZZLE DESIGN
grain.ae = grain.aeat * grain.at; % [mm^2]
grain.Lconv = (grain.Re0*2 - sqrt(4*grain.at/pi))/(2*tan(deg2rad(grain.nozzle_conv_angle_degs)));
grain.Ldiv  = (sqrt(4*grain.at*grain.aeat/pi) - sqrt(4*grain.at/pi))/(2*tan(deg2rad(grain.nozzle_div_angle_degs)));
grain.Lnozzle = grain.Lconv + grain.Ldiv;

%% FULL ENGINE PLOT
if plots
%plot_engine_xsection(grain.L0*1e3, grain.Ri0*1e3, grain.Re0*1e3, grain.at*1e6, grain.ae*1e6, grain.Lconv*1e3, grain.Ldiv*1e3, grain.nozzle_conv_angle_degs, grain.nozzle_div_angle_degs);
end

