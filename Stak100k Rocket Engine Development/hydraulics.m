%--------------------------------------------------------------------------
% INITIALIZATION
%--------------------------------------------------------------------------
warning('off', 'all')
restoredefaultpath
PathsInit
load('grain.mat')

%--------------------------------------------------------------------------
% INPUT PARAMETERS
%--------------------------------------------------------------------------

Pend = 20e5;        %[Pa]

rhoInitial  = 998;       %[kg/m^3] 18*C water temperature density
muInitial   = fluid.mu_vec(1);


d1 = 0.05;         %[m] big tube diameter
l_tubes1 = 5;      %[m] Total length of tube sections
tburn  = 25 + 0.1*25;      %[s], + 4secs for safety windows is hardcoded later
BDratio = 3;       %[-]

% Material properties (Stainless Steel 316L)
STank   = 170e6;   %[Pa]
rhoTank = 8000;    %[kg/m^3]

CoolingJacket.NChannels = N_channels;                 % Number of channels
CoolingJacket.BendAngle = grain.nozzle.thetai + grain.nozzle_conv_angle_degs*pi/180;   %[rad] Bend angle of water path!

% Sample distributions along channel length (placeholder values)
CoolingJacket.ReChannels      = fluid.Re;
CoolingJacket.muBulkChannels  = fluid.mu_vec;
CoolingJacket.muWallChannels  = fluid.muWall;
CoolingJacket.rhoChannels     = rhoInitial.*ones(1,max(size(fluid.mu_vec)));
CoolingJacket.z               = fluid.Ds;               % ascissa curvilinea
CoolingJacket.dChannels       = fluid.HydraulicDiameter;

ValvedP=0.1*1e5;
FlowmeterdP = 0; % negligible

Lequivalent = (sqrt(2*grain.at/pi))*2*pi;  % nozzle circumference @ Ax/At = 2 (subsonic location)          

[Pstart,Tank] = IdraulicSizing(Pend,mdot,tburn,BDratio,rhoInitial,muInitial,d1,l_tubes1,STank,rhoTank,...
    CoolingJacket,ValvedP,FlowmeterdP,Lequivalent);

fprintf('\n--- PRESSURANT HYDRAULIC SIZING RESULTS ---\n');
fprintf('Ptank_min = %.2f bar\n', Pstart/1e5);
fprintf('Tank Volume = %.4f m³ | Nitrogen Vol = %.4f m³ |  Mass = %.2f kg\n', ...
    Tank.Volume, Tank.Nvolume, Tank.Mass);
fprintf('-------------------------------\n\n');