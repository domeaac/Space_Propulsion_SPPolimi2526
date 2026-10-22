%% COOLING JACKET TEST FILE 
clear;
clc; 
close all; 


%% init solver 
% nozzle test 
nozzleTest
y_flipped_explicit = @(x) y_Line(-x) .* (x <= -x0 & x > -x1) + ...
                          y_C1(-x)   .* (x <= -x1 & x > -x2) + ...
                          y_C2(-x)   .* (x <= -x2 & x > -x3) + ...
                          y_P(-x)    .* (x <= -x3 & x >= -x4) + Rt;
domain = -[-0.3 0.13];
xValues = linspace(domain(1), domain(2), 10000); 
yValues = y_flipped_explicit(xValues); 
[~ , minIdx] = min(abs(yValues)); 
x_t = xValues(minIdx); 


%%
                        

% constant params for the solver
Const = wallThermo.init.initConst("TotalTemperature", 3000, ... 
                                  "NozzleCurve", y_flipped_explicit,...
                                  "x_throat", x_t,...
                                  "Rc_t", 0.8);
% gas related quantities from CEA
Gas = wallThermo.init.initGas("P0", 70*1e5,...
                              "cstar", 1426,...
                              "Prandtl", 0.6,...
                              "gamma", 1.238,...
                              "mu0", 0.781*1e-4,...
                              "cp", 1.9*1e3);

% Hot wall side imposed 
BC = wallThermo.init.initBC( "T_hotWall",...
                              2200); 

% generate mesh 
nLayers = 2; 
thickness = [15 0.5] * 1e-3; %layer grows from the wall to hot gas
domain = -[-0.3 0.13]; 
k_per_layer = [200 0.9]; 

Mesh = wallThermo.NozzleMesher.mesh(Const.NozzleCurve, domain, nLayers, thickness, 0.02); 
wallThermo.plotWall(Mesh, "InflationFactor", 3); 
Sol = wallThermo.solverHotSide(Mesh, k_per_layer,1,BC,Const,Gas); 


wallThermo.plotSolution(Sol, "thicknessScale", 2)









