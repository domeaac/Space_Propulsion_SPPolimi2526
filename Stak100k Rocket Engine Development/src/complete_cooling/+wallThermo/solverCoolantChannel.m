function fluid = solverCoolantChannel(q_vec, r_mat,NChannels, mdot, T_i, T_panel_vec, x0_i, options)
% NEWNEWCOOLCHANNEL Computes the thermal and fluid properties of a cooling channel.
%
% Inputs:
%   q_vec       [1 x N]   Surface heat flux for each segment [W/m^2]
%   r_mat       [2 x N+1] Nodal coordinates (x, y) of the channel [m]
%                         r_mat(1,:) = axial coordinates x
%                         r_mat(2,:) = radial coordinates y
%   mdot        [1 x 1]   Fluid mass flow rate [kg/s]
%   T_i         [1 x 1]   Fluid inlet temperature [K]
%   T_panel_vec [1 x N]   Panel wall temperature at each segment [K]
%   x0_i        [2 x 1]   Initial guess for solver [velocity; diameter]
%   options     [struct]  Name-value pairs (e.g., 'correlation')
%
% Outputs:
%   fluid       [struct]  Fluid thermophysical properties evaluated per segment

arguments
    q_vec (1,:) double
    r_mat (2,:) double
    NChannels (1,1) double {mustBeInteger}
    mdot (1,1) double
    T_i (1,1) double
    T_panel_vec (1,:) double
    x0_i (3,1) double
    options.correlation (1,1) string {mustBeMember(options.correlation, ["Gnielinski", "DittusBolt"])} = "Gnielinski"; 
end 

rho = 997; 
Mm_water = 0.018015;
N = length(q_vec);

% here there is already appendix for thoose formulations (from mu to dittus bolter and gnielinki)
mu = @(T) 1.856e-14 .* exp(4209 ./ T + 0.04527 .* T - 3.376e-5 .* T.^2);
eta = @(T) -1.48445 + 4.12292*(T / 298.15) - 1.63866*(T / 298.15).^2;
K = @(T) eta(T) * 0.6065;
cp = @(T) ( -203.6060 + 1523.290*(T/1000) - 3196.413*(T/1000).^2 + ...
            2474.455*(T/1000).^3 + 3.855326./((T/1000).^2) ) / Mm_water;

Re = @(T, v, D) (rho .* v .* D) ./ mu(T);
Pr = @(T) (cp(T) .* mu(T)) ./ K(T);
f_fric = @(Re_val) (1.58 * log(Re_val) - 3.28).^(-2);

Nu_Gn = @(Re_val, Pr_val) ((f_fric(Re_val)/8) .* (Re_val - 1000) .* Pr_val) ./ ...
                          (1 + 12.7 .* sqrt(f_fric(Re_val)/8) .* (Pr_val.^(2/3) - 1));
Nu_Db = @(Re_val, Pr_val) 0.0243 .* Re_val.^(0.8) .* Pr_val.^(0.4); 

if options.correlation == "Gnielinski"
    Nu_correlation = Nu_Gn;
    warning("solverCoolantChannel:payAttention", "INconsistency with nick model, here no roughness is assumed")
else
    Nu_correlation = Nu_Db;
end

Nu = @(h, D, T) h .* D ./ K(T);
% Assuming rectangular section Dh = 2ab / (a+b) , a is imposed by Nchannels , b must be computed 
Dh = @(a,b) 2*a*b/(a+b);

% mDot_fun = @(D, v) rho * (D^2) * v;
mDot_fun = @(a,b,v) a*b*rho*v; % correct for rectuangualr 


T_vec = zeros(1, N);
h_vec = zeros(1, N);
D_vec = zeros(1, N);
v_vec = zeros(1, N);

% channel data
midPoints = zeros(2, N);
channelBase = zeros(1, N);
channelHeight = zeros(1, N);
massFlowRateError = zeros(1, N);




T_inlet = T_i; 
y_guess = log(x0_i(:));
opts = optimoptions('fsolve', 'Display', 'off');

for i = 1:N
    x1 = r_mat(1, i); x2 = r_mat(1, i+1);
    y1 = r_mat(2, i); y2 = r_mat(2, i+1);
    T_panel = T_panel_vec(i);
    
    L = sqrt((x2 - x1)^2 + (y2 - y1)^2);
    S = pi * (y1 + y2) * L; 
    
    % compute local channel base len 
    midYNozzle = (y1+y2)/2;
    midXNozzle = (x1+x2)/2;
    midPoints(:,i)= [midXNozzle; midYNozzle];
    localBaseLen = computeLocalBaseLen(midYNozzle, NChannels); 

    % Energy balance to find exit temperature of current segment
    T_eval = T_inlet + (S * q_vec(i)) / (mdot * cp(T_inlet));
    T_vec(i) = T_eval; 
    
    % Simplified convective heat transfer coefficient calculation 
    h_local = q_vec(i) / (T_panel - T_eval); 
    h_vec(i) = h_local; 
    
    res_fun = @(y) compute_transformed_residuals(y, mdot, h_local, T_eval, mDot_fun, Nu, Nu_correlation, Re, Pr, localBaseLen,Dh,NChannels);
    
    [sol_y, ~, exitflag] = fsolve(res_fun, y_guess, opts);
    
    if exitflag <= 0
        error('Solver failed to converge at segment %d. Exit flag: %d', i, exitflag);
    end
    
    sol_x = exp(sol_y);
    v_vec(i) = sol_x(1);
    D_vec(i) = sol_x(2);

    channelHeight(i) = sol_x(3); 
    channelBase(i) = localBaseLen; 

    % compute massflowrate error 
    massFlowRateError(i) = mdot - NChannels * rho * channelHeight(i)*v_vec(i)*channelBase(i); 
    
    y_guess = sol_y; 
    T_inlet = T_eval; % Update inlet temperature for the next segment
end

% struct assingment
fluid.mu_vec                 = mu(T_vec);
fluid.muWall                 = mu(T_panel_vec); 
fluid.K_vec                  = K(T_vec);
fluid.cp_vec                 = cp(T_vec);
fluid.Re                     = Re(T_vec, v_vec, D_vec);
fluid.Pr                     = Pr(T_vec);
fluid.Tprofile               = T_vec; 
fluid.HcoefConv              = h_vec; 
fluid.HydraulicDiameter      = D_vec; 
fluid.Velocity               = v_vec; 
fluid.ChannelBase            = channelBase; 
fluid.ChannelHeight          = channelHeight; 
fluid.MidPoints              = midPoints; 
fluid.massFlowRateError      = massFlowRateError; 
fluid.BaseOverHeight         = channelBase./channelHeight; 
[fluid.Ds, ~]                = computeDs(r_mat); 



end



%% HELPER FUNCTION SECTION 

function F = compute_transformed_residuals(y, mdot, h_local, T_eval, mDot_fun, Nu, Nu_correlation, Re, Pr,localBaseLen,Dh_fun,NChannels)
% COMPUTE_TRANSFORMED_RESIDUALS Evaluates the normalized residual vector using mapped variables.
%
% Inputs:
%   y              - 2x1 array, unconstrained computational state vector [y1; y2] where x = exp(y).
%   mdot           - [1x1] Fluid mass flow rate
%   h_local        - [1x1] Local heat transfer coefficient
%   T_eval         - [1x1] Fluid evaluation temperature
%   mDot_fun       - Function handle for mass flow rate
%   Nu             - Function handle for analytical Nusselt number
%   Nu_correlation - Function handle for empirical Nusselt number
%   Re             - Function handle for Reynolds number
%   Pr             - Function handle for Prandtl number
%
% Outputs:
%   F              - 2x1 array, evaluated non-dimensional residual vector F(exp(y)).

    v = exp(y(1));
    D = exp(y(2));
    b = exp(y(3)); 
    
    Nu_target = Nu_correlation(Re(T_eval, v, D), Pr(T_eval));
    Nu_current = Nu(h_local, D, T_eval);
    
    F = [
        1 - (mDot_fun(localBaseLen,b, v)*NChannels) / mdot;
        1 - Nu_current / Nu_target;
        1 - Dh_fun(localBaseLen, b) / D
    ];
    
end

% FUNCTION TO COMPUTE THE LENGHT OF THE SQUARE CHANNELS 

function baseLen = computeLocalBaseLen(yNozzle, N)

    % compute initial cirmumference 
    localCirc = 2*pi*yNozzle;
    baseLen = localCirc/N; 
end 


function [dS, S] = computeDs(pointMatrix)
% COMPUTEDS Calculates the differential and cumulative curvilinear abscissa.
% Syntax:
%   [dS, S] = computeDs(pointMatrix)
%
% Inputs:
%   pointMatrix - [2 x N] numeric array representing N discrete position 
%                 vectors in a 2D Cartesian frame.
%
% Outputs:
%   dS - [1 x (N-1)] numeric array of the Euclidean norm (chordal distance) 
%        between consecutive points.
%   S  - [1 x N] numeric array of the cumulative curvilinear abscissa 
%        starting from 0 at the first point.

    deltaP = diff(pointMatrix, 1, 2);
    dS = vecnorm(deltaP, 2, 1);
    S = [0, cumsum(dS)];
end

 