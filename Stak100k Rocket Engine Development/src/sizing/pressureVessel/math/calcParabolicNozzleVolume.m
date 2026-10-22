function V = calcParabolicNozzleVolume(a, b, Rt, L)
%==========================================================================
% FUNCTION: calcParabolicNozzleVolume
%
% DESCRIPTION:
%   Computes the exact analytical volume of a solid of revolution generated 
%   by a parabolic profile:
%
%       r(z) = a*z^2 + b*z + Rt
%
%   The volume is obtained by integrating:
%
%       V = ∫ pi * r(z)^2 dz   from z = 0 to z = L
%
%   The coordinate system is defined as:
%     - z = 0 at the throat plane
%     - z-axis aligned with the nozzle axis
%     - r(z) is the local radius
%
% INPUTS:
%   a  - Quadratic coefficient of the parabolic profile [1/m]
%   b  - Linear coefficient of the parabolic profile [-]
%   Rt - Throat radius at z = 0 [m]
%   L  - Axial length of the divergent section [m]
%
% OUTPUTS:
%   V  - Total internal volume of the nozzle [m^3]
%
% NOTES:
%   - Closed-form analytical solution (no numerical integration)
%   - Assumes axisymmetric geometry
%   - Valid for Rao-type bell nozzle approximations
%
%==========================================================================

    % Integral terms from expansion of pi * r(z)^2
    term1 = (a^2 * L^5) / 5;
    term2 = (a * b * L^4) / 2;
    term3 = ((b^2 + 2 * a * Rt) * L^3) / 3;
    term4 = b * Rt * L^2;
    term5 = Rt^2 * L;
    
    % Final volume
    V = pi * (term1 + term2 + term3 + term4 + term5);
    
end