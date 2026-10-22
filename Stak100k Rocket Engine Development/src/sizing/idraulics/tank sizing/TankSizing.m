function [Vtank, Vnitrogen, TotWallVolume, TankMass] = TankSizing (mdot, tburn, BDratio, P, S, rho)
%==========================================================================
% FUNCTION: TankSizing
%
% DESCRIPTION:
%   Computes the required tank size and structural mass for a pressurized
%   propellant system using a blowdown configuration.
%
%   The function evaluates:
%     - Total tank volume
%     - Initial pressurant (nitrogen) volume
%     - Structural wall volume
%     - Tank mass
%
%   The model assumes adiabatic expansion of the pressurant gas and uses
%   simplified structural formulas for preliminary sizing.
%
% INPUTS:
%   mdot    - Liquid mass flow rate [kg/s]
%   tburn   - Burn time [s]
%   BDratio - Blowdown ratio (P_initial / P_final) [-]
%   P       - Internal tank pressure [Pa]
%   S       - Allowable material stress [Pa]
%   rho     - Material density [kg/m^3]
%
% OUTPUTS:
%   Vtank         - Total tank volume [m^3]
%   Vnitrogen     - Initial nitrogen volume [m^3]
%   TotWallVolume - Total structural wall volume [m^3]
%   TankMass      - Tank mass [kg]
%
% ASSUMPTIONS:
%   - Incompressible liquid (constant density)
%   - Adiabatic expansion of pressurant gas
%   - 20% propellant margin for residuals
%   - Cylindrical tank with hemispherical end caps
%   - Fixed aspect ratio: L/D = 3
%   - Thin-wall pressure vessel approximation
%
% NOTES:
%   - External pressure is assumed to be 1 bar
%   - Structural formulas are simplified (ASME-like)
%   - Suitable for preliminary sizing only
%
%==========================================================================


%--------------------------------------------------------------------------
% PRESSURE CORRECTION
%--------------------------------------------------------------------------
% Remove external atmospheric pressure (1 bar)
P = P - 1e5;


%--------------------------------------------------------------------------
% MATERIAL SAFETY FACTOR
%--------------------------------------------------------------------------
% Apply safety factor to allowable stress
S = S / 1.25;


%--------------------------------------------------------------------------
% FLUID PROPERTIES
%--------------------------------------------------------------------------
rhow = 1000;     % Water density [kg/m^3]
gamma = 1.4;     % Nitrogen specific heat ratio [-]


%--------------------------------------------------------------------------
% PROPELLANT VOLUME
%--------------------------------------------------------------------------
% Total required liquid mass (including residual margin)
Mw = 1.2 * mdot * (tburn + 4);

% Convert mass to volume
Vw = Mw / rhow;


%--------------------------------------------------------------------------
% PRESSURANT SIZING (ADIABATIC BLOWDOWN)
%--------------------------------------------------------------------------
% Nitrogen volume from adiabatic relation
Vnitrogen = Vw / (BDratio^(1 / gamma) - 1);

% Total tank volume
Vtank = Vw + Vnitrogen;


%--------------------------------------------------------------------------
% GEOMETRIC SIZING
%--------------------------------------------------------------------------
% Tank composed of:
%   - Cylindrical section
%   - Two hemispherical heads
%
% Assumption: L/D = 3  →  L = 6R
R = (3 * Vtank / (22 * pi))^(1/3);
L = 6 * R;


%--------------------------------------------------------------------------
% STRUCTURAL SIZING
%--------------------------------------------------------------------------
% Thickness of hemispherical heads
t_heads = P * R / (2 * S - 0.2 * P);

% Thickness of cylindrical shell
t_shell = P * R / (S - 0.6 * P);


%--------------------------------------------------------------------------
% WALL VOLUME CALCULATION
%--------------------------------------------------------------------------
% Cylindrical wall volume
CylWallVolume = pi * ((R + t_shell)^2 - R^2) * L;

% Hemispherical heads (equivalent to one full sphere)
SphereWallVolume =  4/3 * pi * ((R + t_heads)^3 - R^3);

% Total structural volume
TotWallVolume = CylWallVolume + SphereWallVolume;


%--------------------------------------------------------------------------
% MASS ESTIMATION
%--------------------------------------------------------------------------
TankMass = TotWallVolume * rho;


%--------------------------------------------------------------------------
% SUMMARY
%--------------------------------------------------------------------------
fprintf('\n========== TANK SIZING SUMMARY ==========\n');
fprintf('--- Water ---\n');
fprintf('  Water mass (w/ margin): %.4f  kg\n',  Mw);
fprintf('  Water volume:           %.6f  m^3\n', Vw);
fprintf('--- Pressurant (N2) ---\n');
fprintf('  Nitrogen volume:             %.6f  m^3\n', Vnitrogen);
fprintf('--- Tank Volume ---\n');
fprintf('  Total tank volume:           %.6f  m^3\n', Vtank);
fprintf('--- Geometry ---\n');
fprintf('  Inner radius R:              %.4f  m\n',   R);
fprintf('  Cylinder length L:           %.4f  m\n',   L);
fprintf('  Shell thickness t_shell:     %.4f  mm\n',  t_shell * 1e3);
fprintf('  Head  thickness t_heads:     %.4f  mm\n',  t_heads * 1e3);
fprintf('--- Structural ---\n');
fprintf('  Cylindrical wall volume:     %.6f  m^3\n', CylWallVolume);
fprintf('  Spherical heads volume:      %.6f  m^3\n', SphereWallVolume);
fprintf('  Total wall volume:           %.6f  m^3\n', TotWallVolume);
fprintf('  Tank mass:                   %.4f  kg\n',  TankMass);
fprintf('=========================================\n\n');

end