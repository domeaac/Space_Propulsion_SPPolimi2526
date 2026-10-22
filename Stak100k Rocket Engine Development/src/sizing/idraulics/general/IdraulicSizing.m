function [Pstart,Tank] = IdraulicSizing(Pend,mdot,tburn,BDratio,rhoInitial,muInitial,d1,l_tubes1,STank,rhoTank,...
    CoolingJacket,ValvedP,FlowmeterdP,Lequivalent)
%--------------------------------------------------------------------------
% IDRAULICSIZING
%--------------------------------------------------------------------------
% Computes the required starting pressure and tank characteristics for a
% propellant feed system, accounting for pressure losses across pipes,
% cooling jacket, bends, and components.
%
% INPUTS:
%   Pend            - Required pressure at injector outlet [Pa]
%   mdot            - Mass flow rate [kg/s]
%   tburn           - Burn time [s]
%   BDratio         - Blowdown ratio [-]
%   rhoInitial      - Fluid density at inlet [kg/m^3]
%   muInitial       - Fluid dynamic viscosity at inlet [Pa·s]
%   d1              - Diameter of main feed tube [m]
%   l_tubes1        - Length of main feed tube [m]
%   STank           - Tank allowable stress [Pa]
%   rhoTank         - Tank material density [kg/m^3]
%   CoolingJacket   - Struct containing cooling jacket parameters:
%                       .NChannels
%                       .dChannels
%                       .BendAngle
%                       .ReChannels
%                       .muBulkChannels
%                       .muWallChannels
%                       .rhoChannels
%                       .z
%   ValvedP         - Pressure drop across valve [Pa]
%   FlowmeterdP     - Pressure drop across flowmeter [Pa]
%   Lequivalent     - Circumference on nozzle at area ration of 2 [m]
%
% OUTPUTS:
%   Pstart          - Required starting pressure [Pa]
%   Tank            - Struct containing:
%                       .Volume
%                       .Nvolume
%                       .WallsVolume
%                       .Mass
%
% NOTES:
%   - Pressure losses include: pipe friction, tank discharge, cooling
%     jacket distributed and bend losses, fittings, and instrumentation.
%   - A safety factor of 1.2 is applied to total pressure losses.
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% FLOW VELOCITY
%--------------------------------------------------------------------------
% Flow velocity in first tube
v1 = (4 * mdot) / (rhoInitial * pi * d1^2);


%--------------------------------------------------------------------------
% REYNOLDS NUMBER
%--------------------------------------------------------------------------
% Reynolds number in the main tube
re1 = reynolds(rhoInitial, muInitial, v1, d1);


%--------------------------------------------------------------------------
% TUBE LOSSES (DARCY FORMULATION)
%--------------------------------------------------------------------------
% Friction-based loss coefficient
K_tube1 = DarcyStraightTubes(re1, d1) * l_tubes1 / d1;

% Pressure drop in the tube
dP_tube1 = Compute_dPGeneral(K_tube1, rhoInitial, v1);

%--------------------------------------------------------------------------
% TANK (OUTLET) LOSS
%--------------------------------------------------------------------------
% Modeled as an orifice with discharge coefficient
Cd = 0.6;

% Pressure loss due to discharge
dP_tank = 0.5 * rhoInitial * (v1 / Cd)^2;


%--------------------------------------------------------------------------
% COOLING JACKET MODEL
%--------------------------------------------------------------------------
% The cooling jacket is modeled as multiple channels with:
%   - Distributed friction losses (via integration)
%   - Additional bend losses

N = CoolingJacket.NChannels;           % Number of channels
dChannels = CoolingJacket.dChannels;   % Channel diameter [m]
theta = CoolingJacket.BendAngle;       % Bend angle [rad]

ReChannels      = CoolingJacket.ReChannels;
muBulkChannels  = CoolingJacket.muBulkChannels;
muWallChannels  = CoolingJacket.muWallChannels;
rhoChannels     = CoolingJacket.rhoChannels;
z               = CoolingJacket.z;

% Flow per channel
mdotChannels = mdot / N;

% Distributed pressure drop (integrated along channel length)
dPCoolingJack = dPTubesHeatFlux( ...
    ReChannels, muBulkChannels, muWallChannels, ...
    rhoChannels, mdotChannels, dChannels, z);

%--------------------------------------------------------------------------
% BEND LOSSES (COOLING JACKET) - MITRE BEND
%--------------------------------------------------------------------------
% Density approximated as mean value along channels
rhoBend = mean(rhoChannels);
dBend = mean(dChannels);

% Velocity inside channels
vBend = 4 * mdotChannels / (rhoBend * pi * dBend^2);

% Empirical bend loss coefficient (Idelchik correlation)
KBend = 0.95 * (sin(theta/2))^2 + 2.05 * (sin(theta/2))^4;

% Pressure drop due to bend
dP_bend = 0.5 * KBend * rhoBend * vBend^2;

%--------------------------------------------------------------------------
% BENDS AFTER COOLING JACKET
%--------------------------------------------------------------------------
% Two 90° bends at outlet
rhoHot = CoolingJacket.rhoChannels(end);

% Velocity after heating
vHot = 4 * mdot / (rhoHot * pi * d1^2);

% Loss coefficient for 90° bend (Idelchik)
K90Bend = 0.15+0.2;

% Total pressure drop for two bends
dP90Bend = 2 * 0.5 * K90Bend * rhoHot * vHot^2;

%--------------------------------------------------------------------------
% COLLECTOR AND DIVIDER LOSSES
%--------------------------------------------------------------------------

% For the collector/distributor the channels diameters assumed constant ad
% equal to the exit diameter of the cooling jacket.
d2 = dChannels(end);
dPColl_div = dPDivider_Collector(N,mdot,d1,d2,rhoInitial,rhoHot,muInitial,Lequivalent);

%--------------------------------------------------------------------------
% TOTAL PRESSURE DROP
%--------------------------------------------------------------------------
% Includes all contributions + fixed filter losses (1e5 Pa)
dP_tot = ( ...
    dP_tube1 + ...
    dP_tank + ...
    1e5 + ... % Filter losses
    dP_bend + ...
    dPCoolingJack + ...
    dPColl_div + ...
    dP90Bend + ...
    ValvedP + ...
    FlowmeterdP ...
);


%--------------------------------------------------------------------------
% REQUIRED STARTING PRESSURE
%--------------------------------------------------------------------------
% Apply safety factor to total pressure losses
Pstart = (Pend + 1.2 * dP_tot);


%--------------------------------------------------------------------------
% TANK SIZING
%--------------------------------------------------------------------------
% Initial tank pressure from blowdown ratio
PinitialTank = Pstart * BDratio;

% Compute tank properties
[Vtank, Vnitrogen, TotWallVolume, TankMass] = ...
    TankSizing(mdot, tburn, BDratio, PinitialTank, STank, rhoTank);

% Store results in struct
Tank.Volume      = Vtank;
Tank.Nvolume     = Vnitrogen;
Tank.WallsVolume = TotWallVolume;
Tank.Mass        = TankMass;


%--------------------------------------------------------------------------
% PRESSURE LOSSES SUMMARY
%--------------------------------------------------------------------------
fprintf('\n========== HYDRAULIC SIZING — PRESSURE LOSSES ==========\n');
fprintf('  Tube friction loss           : %8.2f  Pa  (%6.2f bar)\n', dP_tube1,      dP_tube1      / 1e5);
fprintf('  Tank discharge loss          : %8.2f  Pa  (%6.2f bar)\n', dP_tank,        dP_tank        / 1e5);
fprintf('  Filter loss (fixed)          : %8.2f  Pa  (%6.2f bar)\n', 1e5,            1.0);
fprintf('  Cooling jacket friction      : %8.2f  Pa  (%6.2f bar)\n', dPCoolingJack,  dPCoolingJack  / 1e5);
fprintf('  Cooling jacket bend loss     : %8.2f  Pa  (%6.2f bar)\n', dP_bend,        dP_bend        / 1e5);
fprintf('  Two 90-deg bends (outlet)    : %8.2f  Pa  (%6.2f bar)\n', dP90Bend,       dP90Bend       / 1e5);
fprintf('  Collector/divider loss       : %8.2f  Pa  (%6.2f bar)\n', dPColl_div,     dPColl_div     / 1e5);
fprintf('  Valve loss                   : %8.2f  Pa  (%6.2f bar)\n', ValvedP,        ValvedP        / 1e5);
fprintf('  Flowmeter loss               : %8.2f  Pa  (%6.2f bar)\n', FlowmeterdP,    FlowmeterdP    / 1e5);
fprintf('  -------------------------------------------------------\n');
fprintf('  Total loss (no SF)           : %8.2f  Pa  (%6.2f bar)\n', dP_tot,         dP_tot         / 1e5);
fprintf('  Total loss (SF = 1.2)        : %8.2f  Pa  (%6.2f bar)\n', 1.2*dP_tot,     1.2*dP_tot     / 1e5);
fprintf('  -------------------------------------------------------\n');
fprintf('  Required starting pressure   : %8.2f  Pa  (%6.2f bar)\n', Pstart,         Pstart         / 1e5);
fprintf('=========================================================\n\n');

end