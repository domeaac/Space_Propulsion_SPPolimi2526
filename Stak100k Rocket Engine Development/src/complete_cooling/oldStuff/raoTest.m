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

% rao mesh 
nLayers = 3; 
layerThickness = [1 1 1] *1e-2; 

MeshRao = wallThermo.meshRao(alpha, beta, At, Ae, Rcc, thickness, N, Lcc, nLayers, layerThickness ,"PlotNozzle", true);

wallThermo.plotWall(MeshRao);

% now that mesh iscomputed , attempt solution 
