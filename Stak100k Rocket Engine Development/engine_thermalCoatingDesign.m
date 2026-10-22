clear; 
warning('off', 'all')
restoredefaultpath

% Load grain for design  
load(fullfile('data', 'grain.mat'));
grain.ae = grain.aeat*grain.at;
addpath(fullfile('src', 'complete_cooling'))
addpath(fullfile('src', 'complete_cooling', 'data'))

%% GENERATE RAO MESH  
Thermal_Ch_thickness = 0; % placeholder 
CombChamberRadius = grain.Re0; 
alpha_div = deg2rad(grain.nozzle_div_angle_degs); 
beta_conv = deg2rad(grain.nozzle_conv_angle_degs); 
Ae = grain.ae; 
At = grain.at; 
N_points = 1000;          % mesh density 
CombChamberLen = 1;      % placeholder 
RaoThickness = 0.01; 
N_layers = 1; 
layerThickness = [1]*1e-3; %#ok<NBRAK1> % mm 

RaoMesh = wallThermo.NozzleMesher.meshRao(alpha_div, beta_conv, At, Ae,...
                                          CombChamberRadius,RaoThickness,...
                                          N_points,CombChamberLen,N_layers,...
                                          layerThickness,"PlotNozzle",false ); 


ARatio_lim_min = 2; 
ARatio_lim_max = 2; 
ARatio_x = pi*RaoMesh.curveNodes(2,:).^2 / RaoMesh.At; 
firstARatio = ARatio_x(1); 
lastARatio = ARatio_x(end); 

cut1 = [firstARatio, ARatio_lim_min];    % sub
cut2 = [ARatio_lim_min, ARatio_lim_max]; % subSup
cut3 = [ARatio_lim_max, lastARatio];     % sup

% cut mesh
coolJacketMesh = wallThermo.NozzleMesher.cropMesh(RaoMesh, cut2, "sub-sup"); 

% get extreme point of the mesh 
x_left = coolJacketMesh.curveNodes(1,1); 
x_right = coolJacketMesh.curveNodes(1,end); 


%% TEST TO CHECK IF TBC IS  ACTUALLY NEEDED 

k_Jacket = 260; % assuming same material as nozzle 
t_Jacket = 2*1e-3; % m

k_Copper = 260; 
t_Copper = 0.01; % m


innerThickness = []; % grows inward (empty , no innerlayers)
outerThickness = [t_Copper t_Jacket]; % grows outward (only nozzle and channels)

k_Layers = [k_Jacket k_Copper ]; % outermost to innermost 

% edit coolJacketMesh 
noTBCMesh = wallThermo.NozzleMesher.newSetLayers(coolJacketMesh , innerThickness,outerThickness); 

% Get cold wall temperature profile 
COLD_WALL_TEMP = 455.15; 
TColdFunc = wallThermo.utilis.getParabolaTwoPoints(x_left , COLD_WALL_TEMP , x_right, COLD_WALL_TEMP, 1, "ZeroDerivativeAtX1", true); 

% call solver (COLD SOLUTION)
SolutionNoTBC = wallThermo.solverColdSideSigmaCorr(noTBCMesh,k_Layers,grain , "BcTColdFuncHandle", TColdFunc);

% plot 
wallThermo.plotWall(noTBCMesh,"ThicknessScales",[1,0.4] , "sectionTitle", "Cooling Jacket Mesh (No TBC)"); 
% wallThermo.plotWall(coolJacketMesh,"ThicknessScales",[1,1,1,1,1] , "sectionTitle", "Cool Jacket Mesh"); 
wallThermo.plotSolution(SolutionNoTBC, "thicknessScales", [1,0.4],"sectionTitle", "Cooling Jakcet (No TBC)");

%% SOLVE FOR COOLING CHANNEL (with NO TBC)

HeatFlux = SolutionNoTBC.HeatFluxProfile; 
T_coldWall = SolutionNoTBC.TemperatureProfile(:,end); 
HeatExchangerNodes = SolutionNoTBC.Mesh.layerNodes(:,:,1); %outer layer 
mdot = 2; 
initialTemp = 291.15; 
N_channels = 180; 

% initial condition is for hydr diamater , velocity anc channel height 
%fluid = wallThermo.solverCoolantChannelFminCon(HeatFlux, HeatExchangerNodes,N_channels, mdot, initialTemp, T_coldWall,[20 , 0.05, 0.01] ,correlation="Gnielinski");
fluidNoTBC = wallThermo.solverCoolantChannel(HeatFlux, HeatExchangerNodes,N_channels, mdot, initialTemp, T_coldWall,[100 , 0.05, 0.01] ,correlation="DittusBolt");

wallThermo.plotChannel( HeatExchangerNodes, fluidNoTBC,"title_prefix", "Fluid (no TBC): ");



%% COOLING JACKET SECTION (With TBC)

% From wall to hot gas : K YSZ 0.9 K_RSZ 0.6 totThick = 350/2 micrometers

k_Jacket = 260; % assuming same material as nozzle 
t_Jacket = 2*1e-3; % m

k_Copper = 260; 
t_Copper = 0.01; % m

k_intermediate = 0.9; 
t_intermediate = 0.24*1e-3 * 1/2; % m


k_topCoat = 0.6;
t_topCoat = t_intermediate; 

k_bondCoat = wallThermo.utilis.BondCoat(770); % good eyeball, it reaches approx 770K  
t_bondCoat = 0.06*1e-3; % m

innerThickness = [t_bondCoat t_intermediate t_topCoat]; % grows inward 
outerThickness = [t_Copper t_Jacket]; % grows outward 

k_Layers = [k_Jacket k_Copper k_bondCoat k_intermediate k_topCoat]; % outermost to innermost 

% edit coolJacketMesh 
coolJacketMesh = wallThermo.NozzleMesher.newSetLayers(coolJacketMesh , innerThickness,outerThickness); 

% Get cold wall temperature profile 
COLD_WALL_TEMP = 455.15; 
TColdFunc = wallThermo.utilis.getParabolaTwoPoints(x_left , COLD_WALL_TEMP , x_right, COLD_WALL_TEMP, 1, "ZeroDerivativeAtX1", true); 

% call solver (COLD SOLUTION)
SolutionCoolJacket = wallThermo.solverColdSideSigmaCorr(coolJacketMesh,k_Layers,grain , "BcTColdFuncHandle", TColdFunc);

% plot 
wallThermo.plotWall(coolJacketMesh,"ThicknessScales",[1,0.4,25,10,10] , "sectionTitle", "Cooling Jacket Mesh (with TBC)"); 
%wallThermo.plotWall(coolJacketMesh,"ThicknessScales",[1,1,1,1,1] , "sectionTitle", "Cool Jacket Mesh"); 
wallThermo.plotSolution(SolutionCoolJacket, "thicknessScales", [1,0.4,25,10,10],"sectionTitle", "Cooling Jakcet (with TBC)");


%% SOLVE FOR COOLING CHANNEL (with TBC)

HeatFlux = SolutionCoolJacket.HeatFluxProfile; 
T_coldWall = SolutionCoolJacket.TemperatureProfile(:,end); 
HeatExchangerNodes = SolutionCoolJacket.Mesh.layerNodes(:,:,1); %outer layer 
mdot = 2; 
initialTemp = 291.15; 
N_channels = 180; 

% initial condition is for hydr diamater , velocity anc channel height 
%fluid = wallThermo.solverCoolantChannelFminCon(HeatFlux, HeatExchangerNodes,N_channels, mdot, initialTemp, T_coldWall,[20 , 0.05, 0.01] ,correlation="Gnielinski");
fluid = wallThermo.solverCoolantChannel(HeatFlux, HeatExchangerNodes,N_channels, mdot, initialTemp, T_coldWall,[100 , 0.05, 0.01] ,correlation="DittusBolt");

wallThermo.plotChannel( HeatExchangerNodes, fluid, "title_prefix", "Fluid (with TBC): ");

%% ABLATIVE k 0.64 0.71 (S12666) rateOfBurn 0.029 0.033 mm/s 

 burnTime = 25; %s 
 rateOfBurn = 0.029; %mm/s 

% % adding virgin material to take perpendicullar diffusivity into account 
 alpha_perp = 0.0036; % cm^2 / s 
 safetFactor = 1.5; 

% % this assuming a wall heat flux on 750 W / cm^2 
 [ablativeThickness, initialThickness, virginThickness]=ablativeMaterial.ComputeAblativeThickness(burnTime, rateOfBurn, alpha_perp, safetFactor); 
