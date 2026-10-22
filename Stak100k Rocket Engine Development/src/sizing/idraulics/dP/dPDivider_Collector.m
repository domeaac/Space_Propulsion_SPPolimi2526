function dP = dPDivider_Collector(N,mdot,d1,d2,rho,rhoHot,mu,Lequivalent)

%==========================================================================
% FUNCTION: dPDivider_Collector
%
% DESCRIPTION:
%   Computes the total pressure drop across a flow divider and collector
%   system, modeled as a ring distributing flow into multiple channels and
%   then recombining it.
%
% INPUTS:
%   N            - Number of channels (must be even) [-]
%   mdot         - Total mass flow rate [kg/s]
%   d1           - Diameter of main pipe [m]
%   d2           - Diameter of individual channels [m]
%   rho          - Fluid density at inlet conditions [kg/m^3]
%   rhoHot       - Fluid density at outlet (heated) conditions [kg/m^3]
%   mu           - Fluid dynamic viscosity [Pa·s]
%   Lequivalent  - Equivalent length of the collector ring [m]
%
% OUTPUT:
%   dP           - Total pressure drop (divider + collector) [Pa]
%
% NOTES:
%   - The system is modeled as a circular manifold (ring collector)
%   - Flow splits evenly into N channels and recombines symmetrically
%   - Empirical correlations (Idelchik handbook) are used for local losses
%
% ASSUMPTIONS:
%   - N must be even (symmetry of the system)
%   - Constant mass flow rate in each channel
%   - Curvature effects neglected (Dean number assumed small)
%   - Turbulent mixing losses during recombination are neglected
%   - Pressure losses occur only in the fluid being diverted/collected
%
%==========================================================================


% Equivalent length of the ring manifold [m]
% N must be even

vInlet = (4 * mdot) / (rho * pi * d1^2);
mdotChannels = mdot / N;
vChannels = (4 * mdotChannels) / (rho * pi * d2^2);
mdotBranch = mdot / 2; % Mass flow rate through each half of the collector
lsection = Lequivalent / N;

% Model:
% Ring collector distributing flow from a larger pipe into N smaller channels
%
% Assumptions:
% - Constant mass flow rate in each channel
% - Curvature effects neglected (large radius → small Dean number)
% - Turbulent effects due to flow recombination are not considered
% - Pressure losses occur only in the fluid being diverted/collected
%
% Pressure losses considered:
% - Friction along main pipe walls
% - Entrance/bend losses into channels (based on upstream velocity)
% - Losses due to flow splitting


%----------------------------- LOSSES IN MAIN PIPE ------------------------%

q = zeros(N/2,1);
velMain = zeros(N/2,1);
Re = zeros(N/2,1);
f = zeros(N/2,1);

q(1) = mdotBranch;
velMain(1) = (4 * mdotBranch) / (rho * pi * d1^2);
Re(1) = reynolds(rho,mu,velMain(1),d1);
f(1) = DarcyStraightTubes(Re(1), d1);

dPSection = 0.5 * rho * velMain(1)^2 * f(1) * lsection / d1;

for i = 1:floor(N/2)-1
    
    q(i+1) = q(i) - mdotChannels;
    
    velMain(i+1) = (4 * q(i+1)) / (rho * pi * d1^2);
    Re(i+1) = reynolds(rho,mu,velMain(i+1),d1);
    f(i+1) = DarcyStraightTubes(Re(i+1), d1); 

    dPSection = dPSection + 0.5 * rho * velMain(i+1)^2 * f(i+1) * lsection / d1;

end


%============================== DIVIDER ===================================

% Tee junction losses at channel entrances
%
% Reference: handbook page 254 (welded tees)
% Main-stream velocity approximated using upstream section velocity
%
% Correlation from chart (page 280), valid for wb/wc > 2

dPChannels = 0;

% for i = 1:floor(N/2)
%     K(i) = 0.9 * (0.34 + (vChannels / velMain(i))^2);
%     dPChannels = dPChannels + 0.5 * K(i) * rho * velMain(i)^2;
% end
K(i) = 0.9 * (0.34 + (vChannels / velMain(i))^2);
dPChannels = dPChannels + 0.5 * K(1) * rho * velMain(1)^2;

%--------------------- LOSSES OF FIRST T DIVIDER --------------------------%

% Curved separator with r0/dc = 2 (handbook page 290)
KSeparator = 0.2;
dPSeparator = 0.5 * KSeparator * rho * vInlet^2;


% Total divider losses
dPDivider = dPSeparator + dPChannels + dPSection;


%============================= COLLECTOR ==================================

% Same geometry as divider but reversed flow direction

% Friction losses scaled for density variation
dPSectionHot = dPSection * rho / rhoHot;


%--------------------- FLOW RECOMBINATION (TEE) ---------------------------%

% Same curved tee model as divider, used as junction
% Equal flow recombination assumed
% Empirical coefficient from interpolation

KCombiner = 0.3;

vOutlet = (4 * mdot) / (rhoHot * pi * d1^2);

dPCombiner = 0.5 * rhoHot * KCombiner * vOutlet^2;


%--------------------- CONVERGING TEE LOSSES ------------------------------%

K = -0.5; 

dPChannelsOut = 0.5 * (K * rhoHot * (velMain(1) * rho / rhoHot)^2);

% Recombined flow assumed fully turbulent (loss independent of Re)


%--------------------------- TOTAL COLLECTOR ------------------------------%

dPCollector = dPSectionHot + dPChannelsOut + dPCombiner;


%----------------------------- TOTAL dP -----------------------------------%

dP = dPDivider + dPCollector;

end