function fluid = solverCoolantChannelFminCon(q_vec, r_mat, NChannels, mdot, T_i, T_panel_vec, x0_i, options)
% SOLVERCOOLANTCHANNELFMINCON Computes thermal and fluid properties of a cooling channel via constrained optimization.
%
% Inputs:
%   q_vec       [1 x N]   Surface heat flux for each segment [W/m^2]
%   r_mat       [2 x N+1] Nodal coordinates (x, y) of the channel [m]
%   NChannels   [1 x 1]   Total number of parallel cooling channels
%   mdot        [1 x 1]   Total fluid mass flow rate [kg/s]
%   T_i         [1 x 1]   Fluid inlet temperature [K]
%   T_panel_vec [1 x N]   Panel wall temperature at each segment [K]
%   x0_i        [3 x 1]   Initial guess [velocity; hydraulic_diameter; channel_height]
%   options     [struct]  Name-value pairs (e.g., 'correlation')
%
% Outputs:
%   fluid       [struct]  Fluid thermophysical properties evaluated per segment

arguments
    q_vec (1,:) double
    r_mat (2,:) double
    NChannels (1,1) double {mustBeInteger, mustBePositive}
    mdot (1,1) double
    T_i (1,1) double
    T_panel_vec (1,:) double
    x0_i (3,1) double
    options.correlation (1,1) string {mustBeMember(options.correlation, ["Gnielinski", "DittusBolt"])} = "Gnielinski"; 
end 

rho = 997; 
Mm_water = 0.018015;
N = length(q_vec);

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
else
    Nu_correlation = Nu_Db;
end

Nu = @(h, D, T) h .* D ./ K(T);
Dh = @(a,b) 2*a*b/(a+b);
mDot_fun = @(a,b,v) a*b*rho*v; 

T_vec = zeros(1, N);
h_vec = zeros(1, N);
D_vec = zeros(1, N);
v_vec = zeros(1, N);

midPoints = zeros(2, N);
channelBase = zeros(1, N);
channelHeight = zeros(1, N);
massFlowRateError = zeros(1, N);

T_inlet = T_i; 
x_guess = x0_i(:);

opts = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point');
lb = [1e-3; 1e-6; 1e-6]; 
ub = [100; 1; 1];       

for i = 1:N
    x1 = r_mat(1, i); x2 = r_mat(1, i+1);
    y1 = r_mat(2, i); y2 = r_mat(2, i+1);
    T_panel = T_panel_vec(i);
    
    L = sqrt((x2 - x1)^2 + (y2 - y1)^2);
    S = pi * (y1 + y2) * L; 
    
    midYNozzle = (y1+y2)/2;
    midXNozzle = (x1+x2)/2;
    midPoints(:,i) = [midXNozzle; midYNozzle];
    
    localBaseLen = computeLocalBaseLen(midYNozzle, NChannels); 

    % === ASPECT RATIO CONSTRAINT ===
    % Enforce localBase / b <= 10  -->  b >= localBase / 10
    lb(3) = localBaseLen / 8;
    
    % Protect the solver: Ensure initial guess is strictly inside the new bound
    if x_guess(3) < lb(3)
        x_guess(3) = lb(3) * 1.05; % Push it slightly above the boundary
    end
    % ===============================

    T_eval = T_inlet + (S * q_vec(i)) / (mdot * cp(T_inlet));
    T_vec(i) = T_eval; 
    
    h_local = q_vec(i) / (T_panel - T_eval); 
    h_vec(i) = h_local; 
    
    obj_fun = @(x) compute_objective(x, mdot, h_local, T_eval, mDot_fun, Nu, Nu_correlation, Re, Pr, localBaseLen, Dh, NChannels);
    
    % Call fmincon with the dynamically updated lb array
    [sol_x, fval, exitflag] = fmincon(obj_fun, x_guess, [], [], [], [], lb, ub, [], opts);
    
    if exitflag <= 0 || fval > 1e-4
        error('Solver failed to converge at segment %d. Exit flag: %d, Objective Value: %.4e', i, exitflag, fval);
    end
    
    v_vec(i)         = sol_x(1);
    D_vec(i)         = sol_x(2);
    channelHeight(i) = sol_x(3); 
    channelBase(i)   = localBaseLen; 

    massFlowRateError(i) = mdot - NChannels * rho * channelHeight(i) * v_vec(i) * channelBase(i); 
    
    x_guess = sol_x; 
    T_inlet = T_eval; 
end

fluid.mu_vec            = mu(T_vec);
fluid.K_vec             = K(T_vec);
fluid.cp_vec            = cp(T_vec);
fluid.Re                = Re(T_vec, v_vec, D_vec);
fluid.Pr                = Pr(T_vec);
fluid.Tprofile          = T_vec; 
fluid.HcoefConv         = h_vec; 
fluid.HydraulicDiameter = D_vec; 
fluid.Velocity          = v_vec; 
fluid.ChannelBase       = channelBase; 
fluid.ChannelHeight     = channelHeight; 
fluid.MidPoints         = midPoints; 
fluid.massFlowRateError = massFlowRateError; 
fluid.BaseOverHeight    = channelBase ./ channelHeight; 

end

%% HELPER FUNCTION SECTION 

function J = compute_objective(x, mdot, h_local, T_eval, mDot_fun, Nu, Nu_correlation, Re, Pr, localBaseLen, Dh_fun, NChannels)
% COMPUTE_OBJECTIVE Calculates the sum of squared non-dimensional residuals for fmincon.
%
% Inputs:
%   x              - [3x1] Physical state vector [v; D; b]
%   mdot           - [1x1] Total fluid mass flow rate [kg/s]
%   h_local        - [1x1] Local heat transfer coefficient [W/m^2K]
%   T_eval         - [1x1] Fluid evaluation temperature [K]
%   mDot_fun       - Function handle for mass flow rate
%   Nu             - Function handle for analytical Nusselt number
%   Nu_correlation - Function handle for empirical Nusselt number
%   Re             - Function handle for Reynolds number
%   Pr             - Function handle for Prandtl number
%   localBaseLen   - [1x1] Computed channel base width (a) [m]
%   Dh_fun         - Function handle for hydraulic diameter
%   NChannels      - [1x1] Total number of cooling channels
%
% Outputs:
%   J              - [1x1] Scalar objective value (L2 norm squared of residuals)

    v = x(1);
    D = x(2);
    b = x(3); 
    
    Nu_target = Nu_correlation(Re(T_eval, v, D), Pr(T_eval));
    Nu_current = Nu(h_local, D, T_eval);
    
    F1 = 1 - (mDot_fun(localBaseLen, b, v) * NChannels) / mdot;
    F2 = 1 - Nu_current / Nu_target;
    F3 = 1 - Dh_fun(localBaseLen, b) / D;
    
    J = F1^2 + F2^2 + F3^2;
end

function baseLen = computeLocalBaseLen(yNozzle, N)
% COMPUTELOCALBASELEN Computes the circumferential width of a single channel.
%
% Inputs:
%   yNozzle - [1x1] Radial coordinate of the nozzle
%   N       - [1x1] Number of channels
%
% Outputs:
%   baseLen - [1x1] Width dimension a

    localCirc = 2 * pi * yNozzle;
    baseLen = localCirc / N; 
end