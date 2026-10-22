clc, clear, close all

% DATA INPUT
% BARIA Bates-motor
% Grain
% Geometric properties
wt = 30;                              % [mm] web thickness
L0 = 290;                             % [mm]
Ri0 = 50;                             % [mm]
Re0 = Ri0 + wt;                       % [mm]
V = pi*(Re0^2 - Ri0^2)*L0;            % [mm^3] initial total grain volume

A_t = pi.*[28.80,25.26,21.81].^2 ./4;       % [mm^2]

% Mass Properties:
mperc_vec = [0.68, 0.18, 0.14];             % [] OX, FU, BI percentages
rho_vec   = [1950, 2700, 920]*1e-9;         % [kg/mm^3]
rho_av = 1/sum(mperc_vec./rho_vec);         % [kg/mm^3] average propellant density
M_tot = rho_av*V;                           % [kg]

% Experimental data loading: 9 batches with 3 motor firings each (@ LP,MP,HP)
data = load("tracesbar1.mat");
batch = fieldnames(data);                   % names of batches fields of the struct

%% PLOTTING THE LOADED DATA

plotLoadedData

%% TIME INTERVALS

timeIntervals

%% Performance of C.C.:
% experimental characteristic velocity:
% Under assumption that c is changing only negligibly wrt the pressure, we
% can average the entire dataset (assuming all variables for all batches to be indipendent and identically
% distributed (same history of production))

% cstar mean value
c = mean(c_vec);
% c = 1.5128*1e3;
% impulse mean value
Is = c/9.80665;

% Stima Deviazione Standard (Deviazione std campionaria (Bessel))
Sigma_c = sqrt(1/(27-1) * sum((c_vec - c).^2)); % [m/s]

%% uncertaintes...
UncertaintyScript
figure;
% plot data pre-regression 
plot(log(p_val),log(r_val),'o'); hold on;
% plot regression outcome
plot([min(log(p_val)),max(log(p_val))], polyval([n,log(a)],[min(log(p_val)),max(log(p_val))]));

%% MOTOR SIMULATION
MotorGiorgio

%% FINAL SUMMARY PLOTS
FinalPlots
