%==========================================================================
% SCRIPT: Pressure Drop Calculation in Flow System
%
% DESCRIPTION:
%   Computes the total pressure drop in a fluid system including:
%     - Straight tube friction losses
%     - Tank discharge loss
%     - Cooling jacket distributed losses
%     - Bend (minor) losses
%
%   The model assumes incompressible flow and constant properties within
%   each section.
%
% INPUT PARAMETERS:
%   Pend          - Final pressure [Pa]
%   mdot          - Mass flow rate [kg/s]
%   rho           - Fluid density [kg/m^3]
%   mu            - Dynamic viscosity [Pa·s]
%   d1            - Main tube diameter [m]
%   l_tubes1      - Total tube length [m]
%   tburn         - Burn time [s]
%   BDratio       - Blowdown ratio [-]
%
% OUTPUT:
%   Pstart        - Required starting pressure [Pa]
%
% NOTES:
%   - Assumes incompressible flow (constant density)
%   - Includes both distributed and lumped losses
%   - Cooling jacket modeled as parallel channels
%   - Some components (valves, filters) are approximated
%
%==========================================================================

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
mdot = 0.1;        %[kg/s]
rho  = 1000;       %[kg/m^3]
mu   = 1e-3;       %[Pa·s]

d1 = 0.01;         %[m]
l_restriction = 0.03; %[m] (currently unused)
l_tubes1 = 1;      %[m] Total length of tube sections

tburn  = 500;      %[s]
BDratio = 3;       %[-]

% Material properties (Stainless Steel 316L)
STank   = 170e6;   %[Pa]
rhoTank = 8000;    %[kg/m^3]

%--------------------------------------------------------------------------
% FLOW VELOCITY
%--------------------------------------------------------------------------
% Continuity equation:
% v = (4 * mdot) / (rho * pi * d^2)
v1 = (4 * mdot) / (rho * pi * d1^2);

%--------------------------------------------------------------------------
% REYNOLDS NUMBER
%--------------------------------------------------------------------------
% Reynolds number in the main tube
re1 = reynolds(rho, mu, v1, d1);

%--------------------------------------------------------------------------
% TUBE LOSSES (DARCY FORMULATION)
%--------------------------------------------------------------------------
% Friction-based loss coefficient
K_tube1 = DarcyStraightTubes(re1, d1) * l_tubes1 / d1;

% Pressure drop in the tube
dP_tube1 = Compute_dPGeneral(K_tube1, rho, v1);

%--------------------------------------------------------------------------
% TANK (OUTLET) LOSS
%--------------------------------------------------------------------------
% Modeled as an orifice with discharge coefficient
Cd = 0.8;

% Pressure loss due to discharge
dP_tank = 0.5 * rho * (v1 / Cd)^2;

%--------------------------------------------------------------------------
% COOLING JACKET MODEL
%--------------------------------------------------------------------------
% The cooling jacket is modeled as multiple channels with:
%   - Distributed friction losses (via integration)
%   - Additional bend losses

N = 10;                 % Number of channels
dChannels = 0.01;       %[m]
LChannels = 3;          %[m] (not used in current model)
theta = deg2rad(120);   %[rad] Bend angle

% Sample distributions along channel length (placeholder values)
ReChannels      = [1; 2; 3];
muBulkChannels  = [1; 2; 3];
muWallChannels  = [1; 2; 3];
rhoChannels     = [1; 2; 3];
z               = [1; 2; 3];

% Flow per channel
mdotChannels = mdot / N;

% Velocity in each channel
vChannels = (4 * mdotChannels) / (rho * pi * dChannels^2);

% Distributed pressure drop (integrated along z)
dPCoolingJack = dPTubesHeatFlux( ...
    ReChannels, muBulkChannels, muWallChannels, ...
    rhoChannels, vChannels, dChannels, z);

%--------------------------------------------------------------------------
% BEND LOSSES (COOLING JACKET)
%--------------------------------------------------------------------------
% Density approximated as mean value along channels
rhoBend = mean(rhoChannels);

% Empirical bend loss coefficient (Idelchik correlation)
KBend = 0.95 * (sin(theta/2))^2 + 2.05 * (sin(theta/2))^4;

% Pressure drop due to bend
dP_bend = 0.5 * KBend * rhoBend * vChannels^2;

%--------------------------------------------------------------------------
% BENDS AFTER COOLING JACKET
%--------------------------------------------------------------------------
% Two 90° bends (mitre approximation)
K90Bend = 0.21;

dP90Bend = 2 * 0.5 * K90Bend * rho * v1^2;

%--------------------------------------------------------------------------
% TOTAL PRESSURE DROP
%--------------------------------------------------------------------------
% Includes estimated additional lumped losses (~1e5 Pa)
dP_tot = ( ...
    dP_tube1 + ...
    dP_tank + ...
    1e5 + ...
    dP_bend + ...
    dPCoolingJack + ...
    dP90Bend ...
);

%--------------------------------------------------------------------------
% REQUIRED STARTING PRESSURE
%--------------------------------------------------------------------------
% Safety factor applied to total losses
Pstart = (Pend + 1.2 * dP_tot);

%--------------------------------------------------------------------------
% TANK SIZING
%--------------------------------------------------------------------------
PinitialTank = Pstart * BDratio;

[Vtank, Vnitrogen, TotWallVolume, TankMass] = ...
    TankSizing(mdot, tburn, BDratio, PinitialTank, STank, rhoTank);