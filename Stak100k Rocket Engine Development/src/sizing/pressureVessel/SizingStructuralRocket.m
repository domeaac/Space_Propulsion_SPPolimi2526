function [m_tot, D_in_esterno, t] = SizingStructuralRocket( ...
    Pc, Pt, nozzleGeometry, Length, S, rho, S_nozzle, rho_nozzle, extra)

%==========================================================================
% FUNCTION: SizingStructuralRocket
%
% DESCRIPTION:
%   Performs structural sizing and mass estimation of a solid rocket motor.
%   The following components are included:
%       - Cylindrical case
%       - Hemispherical head
%       - Convergent-divergent nozzle (Rao bell approximation)
%
%   The nozzle thickness is computed assuming conical shell behavior under
%   pressure, while its volume is evaluated using a Rao parabolic profile.
%   This leads to a conservative (overestimated) nozzle mass.
%
% INPUTS:
%   Pc              - Chamber pressure [Pa]
%   Pt              - Throat pressure [Pa]
%   nozzleGeometry  - Geometry vector:
%                       [D_in, D_out, Dt, alpha_conv, alpha_div, a, b, RaoL]
%                       D_in       - Internal case diameter [m]
%                       D_out      - Nozzle exit diameter [m]
%                       Dt         - Throat diameter [m]
%                       alpha_conv - Convergent half-angle [rad]
%                       alpha_div  - Divergent half-angle [rad]
%                       a, b       - Rao profile parameters [-]
%                       RaoL       - Nozzle length [m]
%   Length          - Cylindrical section length [m]
%   S               - Allowable stress (case material) [Pa]
%   rho             - Density of case material [kg/m^3]
%   S_nozzle        - Allowable stress (nozzle material) [Pa]
%   rho_nozzle      - Density of nozzle material [kg/m^3]
%
% OUTPUTS:
%   m_tot           - Total structural mass [kg]
%   D_in_esterno    - External diameter of case [m]
%   t               - Thickness vector [mm]:
%                       t(1) = case thickness
%                       t(2) = head thickness
%                       t(3) = divergent thickness
%                       t(4) = convergent thickness (normal direction)
%
% ASSUMPTIONS:
%   - Safety factor Fs = 1.25 applied to allowable stresses
%   - External pressure = 1 bar (subtracted from Pc and Pt)
%   - Thin-wall pressure vessel approximations
%   - Materials are isotropic and temperature-independent
%   - Perfect welds (efficiency = 1)
%   - Thermal and corrosion effects neglected
%
% NOTES:
%   - Nozzle volume is computed using Rao parabolic profile
%   - Thickness is computed assuming conical shell approximation
%   - Thickness output is converted to mm at the end
%   - Ensure density units are consistent with volume (SI recommended)
%
%==========================================================================


%--------------------------------------------------------------------------
% INITIALIZATION
%--------------------------------------------------------------------------
t  = zeros(3,1);   % Thickness vector: [case, head, divergent]
Fs = 1.15;         % Safety factor

% Extract geometry parameters
D_in       = nozzleGeometry(1);
D_out      = nozzleGeometry(2);
Dt         = nozzleGeometry(3);
alpha_conv = nozzleGeometry(4);
alpha_div  = nozzleGeometry(5);
a          = nozzleGeometry(6);
b          = nozzleGeometry(7);
RaoL       = nozzleGeometry(8);


%--------------------------------------------------------------------------
% PRESSURE CORRECTION
%--------------------------------------------------------------------------
% Convert absolute pressure to gauge pressure (subtract 1 bar)
Pc = Pc - 1e5;
Pt = Pt - 1e5;


%--------------------------------------------------------------------------
% STRUCTURAL SIZING
%--------------------------------------------------------------------------
% Apply safety factor
S         = S / Fs;
S_nozzle  = S_nozzle / Fs;

% Case thickness (thin cylinder approximation)
t(1) = (Pc * D_in / 2) / (S - 0.6 * Pc);

% External diameter of case
D_in_esterno = D_in + 2 * t(1);

% Head thickness (hemispherical shell)
t(2) = (Pc * D_in_esterno / 2) / (2 * S + 0.8 * Pc);

% Divergent thickness (conical shell)
t(3) = (Pt * D_out / 2) / (cos(alpha_div) * (S_nozzle - 0.6 * Pt));

% External throat diameter
Dt_esterno = Dt + 2 * t(3)/cos(alpha_div);


%--------------------------------------------------------------------------
% NOZZLE VOLUME 
%--------------------------------------------------------------------------
LDivergent = (D_out - Dt) / (2 * tan(alpha_div));

VExternal = calcola_volume_tronco_cono(Dt_esterno, rad2deg(alpha_div), LDivergent);
VInternal = calcola_volume_tronco_cono(Dt, rad2deg(alpha_div), LDivergent);

% Structural volume of divergent section
V_div = VExternal-VInternal;
%--------------------------------------------------------------------------
% CYLINDRICAL CASE VOLUME
%--------------------------------------------------------------------------
V_cilindro = pi/4 * (D_in_esterno^2 - D_in^2) * Length;

%--------------------------------------------------------------------------
% HEAD VOLUME (HEMISPHERICAL)
%--------------------------------------------------------------------------
DInternaHead = D_in_esterno - 2 * t(2);

VExternalHead = (2/3) * pi * (D_in_esterno/2)^3;
VInternalHead = (2/3) * pi * (DInternaHead/2)^3;

% Net structural volume
V_testa = VExternalHead - VInternalHead;


%--------------------------------------------------------------------------
% MASS: CASE + HEAD
%--------------------------------------------------------------------------
m_case = rho * (V_testa + V_cilindro);


%--------------------------------------------------------------------------
% CONVERGENT SECTION
%--------------------------------------------------------------------------
% Thickness (conical shell with inclination correction)

if extra
   ddd = sqrt(4*(((pi/4)*Dt^2)*2)/pi);
   t_conv = (Pc * ddd / 2) / ...
         (cos(alpha_conv) * (S_nozzle - 0.6 * Pc));
else
t_conv = (Pc * D_in_esterno / 2) / ...
         (cos(alpha_conv) * (S_nozzle - 0.6 * Pc) +  (Pc / cos(alpha_conv)));
end

% Convert to radial thickness
t_rad_conv = t_conv / cos(alpha_conv);

% External diameters
D_conv_esterno  = D_in_esterno;
Dt_esterno_conv = Dt + 2 * t_rad_conv;

% Internal diameters
D_conv_interno = D_conv_esterno - 2 * t_rad_conv;
Dt_interno     = Dt_esterno_conv - 2 * t_rad_conv;

% Convergent length
L_conv = (D_conv_interno - Dt) / (2 * tan(alpha_conv));

% External volume (frustum)
V_conv_ext = (pi * L_conv / 12) * ...
    (D_conv_esterno^2 + D_conv_esterno * Dt_esterno_conv + Dt_esterno_conv^2);

% Internal volume (frustum)
V_conv_int = (pi * L_conv / 12) * ...
    (D_conv_interno^2 + D_conv_interno * Dt_interno + Dt_interno^2);

% Net structural volume
V_conv = V_conv_ext - V_conv_int;

% Mass of nozzle (convergent + divergent)
m_nozzle = rho_nozzle * ( V_conv + V_div );

% Append convergent thickness
t = [t; t_conv];


%--------------------------------------------------------------------------
% TOTAL MASS
%--------------------------------------------------------------------------
m_tot = m_case + m_nozzle;

% Convert thickness to mm
t = t * 1000;

end