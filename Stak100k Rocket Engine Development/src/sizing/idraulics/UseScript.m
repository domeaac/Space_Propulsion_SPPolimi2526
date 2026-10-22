clear
clc
close all

%--------------------------------------------------------------------------
% INITIALIZATION
%--------------------------------------------------------------------------
% Add necessary paths (user-defined function)
SetPath();

%--------------------------------------------------------------------------
% INPUT PARAMETERS
%--------------------------------------------------------------------------
Pend = 1e5;        %[Pa]
mdot = 1.5;        %[kg/s]
rhoInitial  = 1000;       %[kg/m^3]
muInitial   = 1e-3;       %[Pa·s]
d1 = 0.1;         %[m]
l_tubes1 = 5;      %[m] Total length of tube sections
tburn  = 25;      %[s]
BDratio = 3;       %[-]

% Material properties (Stainless Steel 316L)
STank   = 170e6;   %[Pa]
rhoTank = 8000;    %[kg/m^3]

CoolingJacket.NChannels = 20;                 % Number of channels
CoolingJacket.BendAngle = deg2rad(70);   %[rad] Bend angle

% Sample distributions along channel length (placeholder values)
CoolingJacket.ReChannels      = [1; 2; 3]*1e6;
CoolingJacket.muBulkChannels  = [1; 2; 3];
CoolingJacket.muWallChannels  = [1; 2; 3];
CoolingJacket.rhoChannels     = [1; 2; 3];
CoolingJacket.z               = [1; 2; 3];
CoolingJacket.dChannels       = [0.01; 0.01; 0.01];
ValvedP=0.5*1e5;
FlowmeterdP=0;

Lequivalent=pi;

[Pstart,Tank] = IdraulicSizing(Pend,mdot,tburn,BDratio,rhoInitial,muInitial,d1,l_tubes1,STank,rhoTank,...
    CoolingJacket,ValvedP,FlowmeterdP,Lequivalent)