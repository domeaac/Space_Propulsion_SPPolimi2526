function re = reynolds(rho,mu,v,l)
%==========================================================================
% FUNCTION: reynolds
%
% DESCRIPTION:
%   This function calculates the Reynolds number for a fluid flow, which is
%   a dimensionless quantity used to characterize the flow regime.
%
% INPUTS:
%   rho - Fluid density
%   mu  - Dynamic viscosity of the fluid
%   v   - Flow velocity
%   l   - Characteristic length (e.g., pipe diameter)
%
% OUTPUT:
%   re  - Reynolds number
%
% NOTES:
%   - The Reynolds number is defined as: Re = (rho * v * l) / mu
%   - It is used to distinguish between laminar and turbulent flow:
%       * Re < 2000   → Laminar flow
%       * Re > 4000   → Turbulent flow
%       * Transitional regime in between
%
%==========================================================================

% Compute Reynolds number
re = rho * v * l / mu;

end