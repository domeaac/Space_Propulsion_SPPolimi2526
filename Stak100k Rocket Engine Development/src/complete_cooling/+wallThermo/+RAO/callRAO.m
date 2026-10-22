clear
close all
clc

alpha = deg2rad(15); % divergent angle (conic nozzle)
beta = deg2rad(45); % convergent angle (conic nozzle)
At = 0.0089; % m^2
Ae = 8*At;
Rcc = 0.4;
thickness = 0.03;

N = 100;
Lcc = 1;

[r_nodes, T, N, Lc_Rao,Ld_Rao] = RAO(alpha, beta, Ae, At, Rcc, thickness, N, Lcc, 0); 


