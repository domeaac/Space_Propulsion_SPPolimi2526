%% MOTOR - Mid Pressure Case
At = A_t(2); % MP case;
% a         % mm/s / bar^n
% rho_av    %[kg/mm^3]
% At        % mm^2
% n         % [-]
% Ri0       % [mm]
% Re0       % [mm]
c_star = c*1000;   %mm/s

% CONSERVATION OF MASS ('Pressure Equation') 
dz = @(Ab) ((1e-5)*(a*1e-3)*(rho_av*1e9)*(c_star*1e-3)*Ab./At).^(n/(1-n)); % dz = p^n
     % rho_av [kg/mm^3] -> [kg/m^3]
     % c_star [mm/s] -> [m/s]
     % p [Pa] -> [bar]
     % 

% VIEILLE LAW 
rb = @(p) a*p.^n;
% GRAIN GEOMETRY 
% Ab = @(z) (a*z<(Re0-Ri0)).*(2*pi*(Ri0 + a*z).*(L0-a*z) + pi*(Re0^2 - (Ri0 + a*z).^2));
 Ab = @(z) (a*z<(Re0-Ri0)).*(2*pi*(Ri0 + a*z).*(L0 - 2*a*z) + 2*pi*(Re0^2 - (Ri0 + a*z).^2));

% I.C. for Ab(t):
Ab0 = 2*pi*Ri0*L0 + 2*pi*(Re0^2 - Ri0^2);

% % I.C. for z(t) (should be 0):
% % Solve for z:
% f = @(z) Ab(z) - Ab0;
% z0 = fzero(f, 0);   % 0 = guess iniziale per z in 0, cioè per l'integrale di p^n da 0 a 0: 0
% % check:  = 0? Ok.

t0 = 0; % [s]
tf = 10; % [s]

% dzt = @(t,z) dz(Ab(z));
% [t,z] = ode45(dzt, [t0 tf], z0);

%% RK4
% RUNGE-KUTTA 4 - Generic Implementation

% Define your ODE: dz/dt = f(t, z)
f = @(t, z) dz(Ab(z));  
% Initial conditions
dt = 1e-3;       % passo temporale (riducilo se diverge)
z0 = 0;          % condizione iniziale

% Preallocazione
t_vec = t0 : dt : tf;
N     = length(t_vec);
z_vec = zeros(1, N);
z_vec(1) = z0;

% RK4 loop
for i = 1 : N-1
    t_i = t_vec(i);
    z_i = z_vec(i);

    k1 = f(t_i,        z_i           );
    k2 = f(t_i + dt/2, z_i + dt/2*k1);
    k3 = f(t_i + dt/2, z_i + dt/2*k2);
    k4 = f(t_i + dt,   z_i + dt*k3  );

    z_vec(i+1) = z_i + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end

% Post-processing: ricava p(t) da z(t)
% dz_num = diff(z_vec) ./ dt;
dz_num = dz(Ab(z_vec));
p_sol  = real(dz_num .^ (1/n));

%% Plotting
figure;
subplot(2,1,1)
plot(t_vec, z_vec);
xlabel('t [s]'); ylabel('z [mm·s^{n-1}]');
title('z(t) - RK4')
grid on;

subplot(2,1,2)
plot(t_vec, p_sol);
xlabel('t [s]'); ylabel('p [bar]');
title('p(t) - RK4')
grid on;

figure;
plot(t_vec, Ab(z_vec), 'r', 'LineWidth', 1.5);
xlabel('t [s]');
ylabel('A_b [mm^2]');
title('Burning Area A_b(t)');
grid on;

disp(['a        = ', num2str(a),        ' mm/s/bar^n'])
disp(['rho_av   = ', num2str(rho_av),   ' kg/mm^3'])
disp(['c_star   = ', num2str(c_star/1000),   ' m/s'])
disp(['Ab0      = ', num2str(Ab0),      ' mm^2'])
disp(['At       = ', num2str(At),       ' mm^2'])
disp(['Ab0/At   = ', num2str(Ab0/At)])

p_test = (dz(Ab(0)))^(1/n);
disp(['p_test   = ', num2str(p_test),   ' bar'])

% z=real(z);
% dz_sol = diff(z) ./ diff(t);
% dz_sol=  real(dz_sol);
% p_sol =  dz_sol.^(1/n);

% %% POLT RESULTS OF MOTOR SIMULATION
% figure;
% plot(t,z);
% title('z(t)')
% 
% figure;
% %plot(t(2:end),dz_sol);
% plot([0;t(1:end-1)],[0;p_sol]);
% title('p(t)')
% 
% z_vec = linspace(0, 30, 100);
% figure;
% plot(z_vec, Ab(z_vec));
% xlabel('z [mm]'); ylabel('Ab [mm^2]');
% title('Ab(z)')

