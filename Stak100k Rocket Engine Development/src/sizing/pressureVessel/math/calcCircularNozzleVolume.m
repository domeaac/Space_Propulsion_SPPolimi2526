function V = calcCircularNozzleVolume(Rt, x1, x2)

%==========================================================================
% FUNCTION: calcCircularNozzleVolume
%
% DESCRIPTION:
%   Computes the exact analytical volume of a solid of revolution generated 
%   by a circular throat profile defined as:
%
%       r(x) = Rt + R - sqrt(R^2 - x^2)
%
%   where R = 0.382 * Rt. This represents a convex profile (throat rounding)
%   commonly used in the convergent-divergent transition of rocket nozzles.
%
%   The volume is obtained by the analytical integration of:
%
%       V = ∫ pi * r(x)^2 dx   from x = x1 to x = x2
%
% INPUTS:
%   Rt - Throat radius at the narrowest section (x = 0) [m]
%   x1 - Lower integration limit along the axial coordinate [m]
%   x2 - Upper integration limit along the axial coordinate [m]
%
% OUTPUTS:
%   V  - Total internal volume of the circular section [m^3]
%
% NOTES:
%   - Closed-form analytical solution (no numerical approximation).
%   - The integration limits must satisfy |x/R| <= 1 to avoid complex values.
%   - Assumes axisymmetric geometry around the x-axis.
%
%==========================================================================

    % Define geometric parameters
    R = 0.382 * Rt;
    K = Rt + R; % Translation constant (distance from axis to circle center)

    % Primitive of the polynomial part: ∫ (K^2 + R^2 - x^2) dx
    f_poly = @(x) (K^2 + R^2) * x - (x.^3 / 3);

    % Primitive of the circular/irrational part: ∫ -2 * K * sqrt(R^2 - x^2) dx
    % Note: The factor 2 from the binomial expansion cancels with the 0.5 
    % from the standard integral of the square root.
    f_circ = @(x) - K * (x .* sqrt(R^2 - x.^2) + R^2 * asin(x/R));

    % Total antiderivative function F(x)
    F = @(x) f_poly(x) + f_circ(x);
    
    % Evaluate definite integral: V = pi * [F(x2) - F(x1)]
    V = pi * (F(x2) - F(x1));

end
