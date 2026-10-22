function dP = Compute_dPGeneral(K,rho,v)
%==========================================================================
% FUNCTION: Compute_dPGeneral
%
% DESCRIPTION:
%   This function computes the pressure drop (dP) associated with a given
%   loss coefficient using the dynamic pressure formulation.
%
% INPUTS:
%   K   - Loss coefficient (dimensionless)
%   rho - Fluid density
%   v   - Flow velocity
%
% OUTPUT:
%   dP  - Pressure drop
%
% NOTES:
%   - Based on the standard relation: dP = (1/2) * K * rho * v^2
%   - Commonly used in fluid dynamics to evaluate localized and distributed losses
%
%==========================================================================

% Compute pressure drop using dynamic pressure formula
dP = 0.5 * K * rho * v^2;

end