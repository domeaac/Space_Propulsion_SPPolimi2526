function printGeometricSummary(grain, t_conv2)
%PRINTGEOMETRICSUMMARY  Console table of engine geometry and design-point data.

% --- Derived quantities ---
Dt   = sqrt(4 * grain.at / pi)           * 1e3;   % throat diameter      [mm]
De   = sqrt(4 * grain.at * grain.aeat / pi) * 1e3; % exit diameter        [mm]
Di   = grain.Ri0 * 2                    * 1e3;   % port diameter        [mm]
Dc   = grain.Re0 * 2                    * 1e3;   % casing inner diam.   [mm]
web  = (grain.Re0 - grain.Ri0)          * 1e3;   % web thickness        [mm]
L0mm = grain.L0                         * 1e3;   % segment length       [mm]

Ab0  = pi * Di*1e-3 * grain.L0 + ...
       2 * (pi/4) * ((grain.Re0*2)^2 - (grain.Ri0*2)^2);  % per segment [m^2]
Kn   = Ab0 / grain.at;

LD_port   = grain.L0 / (grain.Ri0 * 2);
LD_casing = grain.L0 / (grain.Re0 * 2);

Ln_mm = grain.nozzle.full_RAO_length * 1e3;   % RAO nozzle length [mm]

t_case  = grain.str.case_thickness                      ;
t_head  = grain.str.head_thickness                      ;
t_divn  = grain.str.divergent_thickness                 ;
t_conv  = grain.str.convergent_thickness_normal_direction ;
De_eng  = grain.str.DeEngine                            ;   % engine OD [mm]
m_str   = grain.str.m_tot;                                        % [kg]

% --- Print ---
fprintf('\n');
fprintf('╔════\n');
fprintf('║          BATES ENGINE – GEOMETRIC SUMMARY            \n');
fprintf('╚════\n\n');
fprintf('║  GRAIN  (per segment, nominal)                       \n');
fprintf('║    Segment length       L0   = %8.2f  mm          ║\n', L0mm);
fprintf('║    Port diameter        Di   = %8.2f  mm          ║\n', Di);
fprintf('║    Casing inner diam.   Dc   = %8.2f  mm          ║\n', Dc);
fprintf('║    Web thickness        w    = %8.2f  mm          ║\n', web);
fprintf('║    L/D  (port)               = %8.3f  [-]         ║\n', LD_port);
fprintf('║    L/D  (casing)             = %8.3f  [-]         ║\n', LD_casing);
fprintf('║    Init. burning area   Ab0  = %8.5f  m²          ║\n', Ab0);
fprintf('╚════\n\n');
fprintf('║  NOZZLE  (RAO parabolic)                             ║\n');
fprintf('║    Throat diameter      Dt   = %8.2f  mm          ║\n', Dt);
fprintf('║    Exit diameter        De   = %8.2f  mm          ║\n', De);
fprintf('║    Area ratio         Ae/At  = %8.3f  [-]         ║\n', grain.aeat);
fprintf('║    RAO nozzle length    Ln   = %8.2f  mm          ║\n', Ln_mm);
fprintf('╚════\n\n');
fprintf('║  STRUCTURE                                           \n');
fprintf('║    Case wall thickness       = %8.2f  mm          \n', t_case);
fprintf('║    Head thickness            = %8.2f  mm          \n', t_head);
fprintf('║    Nozzle div. thickness     = %8.2f  mm          \n', t_divn);
fprintf('║    Nozzle conv. thickness    = %8.2f  mm  (chamber proximity) \n', t_conv);
fprintf('║    Nozzle conv. thickness    = %8.2f  mm  (2 Ae/At  subsonic proximity) \n', t_conv2);
fprintf('║    Engine outer diameter     = %8.2f  m           \n', De_eng);
fprintf('║    Structural mass           = %8.3f  kg          \n', m_str);
fprintf('╚════\n\n');




end
