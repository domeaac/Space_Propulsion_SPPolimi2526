

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