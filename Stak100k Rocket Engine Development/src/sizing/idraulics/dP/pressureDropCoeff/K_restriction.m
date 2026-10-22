function K = K_restriction(d1,d2,l)
%==========================================================================
% FUNCTION: K_restriction
%
% DESCRIPTION:
%   This function calculates the loss coefficient (K) for a flow restriction
%   in a pipe based on the inlet diameter, outlet diameter, and junction length.
%
% INPUTS:
%   d1 - Inlet diameter
%   d2 - Outlet diameter (must be smaller than d1)
%   l  - Length of the restriction (junction length)
%
% OUTPUT:
%   K  - Loss coefficient associated with the restriction
%
% NOTES:
%   - The function assumes a gradual contraction geometry.
%   - Two different empirical formulas are used depending on the angle alpha.
%
%==========================================================================

% d2 exit diameter, d1 inlet diameter, l junction length

% Check that outlet diameter is smaller than inlet diameter
if d2 > d1
    error("oulet must be smaller than inlet");
end

% Calculate the cone angle (alpha)
% Geometry-based relation for contraction angle
alpha = 2 * atan((d1 - d2) / (2 * l));

% Calculate beta parameter (diameter ratio)
beta = d2 / d1;

% Select appropriate empirical formula based on angle
if rad2deg(alpha) <= 45
        % Gradual transition: K depends on sin of half-angle
        K = 0.8 * sin(alpha / 2) * (1 - beta^2);
    elseif rad2deg(alpha) <= 180
        % Steeper transition: Transitioning toward sudden contraction behavior
        % K is influenced by the Borda-Carnot limit (0.5 for alpha = 180)
        K = 0.5 * sqrt(sin(alpha / 2)) * (1 - beta^2)^2;
    else
        % Fallback to sudden contraction for non-physical angles
        K = 0.5 * (1 - beta^2);
end

end