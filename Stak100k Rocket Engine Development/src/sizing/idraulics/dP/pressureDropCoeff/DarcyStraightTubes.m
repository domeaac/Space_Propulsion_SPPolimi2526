function f = DarcyStraightTubes (Re, D)
%==========================================================================
% FUNCTION: DarcyStraightTubes
%
% DESCRIPTION:
%   Computes the Darcy friction factor for flow inside a straight pipe
%   as a function of the Reynolds number.
%
% INPUT:
%   Re - Reynolds number [-]
%
% OUTPUT:
%   f  - Darcy friction factor [-]
%
% METHODOLOGY:
%   The function evaluates the flow regime and applies:
%
%   1) Laminar flow (Re < 2000):
%        Hagen–Poiseuille relation
%        f = 64 / Re
%
%   2) Transitional / moderate turbulent flow (2000 ≤ Re < 1e5):
%        Blasius correlation
%        f = 0.3164 / Re^(0.25)
%
%   3) Fully turbulent flow (Re ≥ 1e5):
%        Haaland correlation (smooth pipe simplification)
%
% ASSUMPTIONS:
%   - Smooth pipe (zero roughness)
%
% NOTES:
%   - The transition region is approximated using the Blasius correlation
%   - Haaland formula is used here without roughness contribution
%
%==========================================================================

%--------------------------------------------------------------------------
% FLOW REGIME IDENTIFICATION
%--------------------------------------------------------------------------

% The Reynolds number determines the applicable correlation
eps = 0.2*1e-3; % in bound value for metal additive manufacturing 
eps_over_d = eps/D;

if Re < 2000

    %----------------------------------------------------------------------
    % LAMINAR FLOW REGIME
    %----------------------------------------------------------------------
    % Hagen–Poiseuille relation:
    % f = 64 / Re
    f = (64 / Re);

% elseif Re < 1e5
% 
%     %----------------------------------------------------------------------
%     % TRANSITIONAL / TURBULENT FLOW REGIME
%     %----------------------------------------------------------------------
%     % Blasius empirical correlation (valid for smooth pipes):
%     % f = 0.3164 / Re^(0.25)
%     f = (0.3164 / Re^(0.25));
% 
else

    %----------------------------------------------------------------------
    % HIGH REYNOLDS TURBULENT REGIME
    %----------------------------------------------------------------------
    % Haaland correlation (simplified for smooth pipes)
    % c is an intermediate variable used in the formulation
    c = -1.8 * log10((6.9 / Re) + (eps_over_d/3.7)^1.11);
    
    % Final Darcy friction factor
    f = (1 / (c^2));

end

end