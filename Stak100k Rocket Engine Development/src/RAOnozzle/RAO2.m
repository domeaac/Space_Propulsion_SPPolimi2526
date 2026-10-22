function [y_patched, r_nodes, tangents, normals, Lc_Rao,Ld_Rao,eta_isoloss, final_theta_i, parabola_params, full_length] = RAO2(alpha, beta, Ae, At, Rcc, thickness, N, Lcc, plots)

 
% Generates the internal wall, the one that is the fluid boundary. The
% thickness is therefore extruded outwards or inwards from this boundary.

% Inputs:
%   alpha      - Divergence semi-angle [rad]
%   beta       - Convergence semi-angle [rad]
%   Ae         - Exit area [m^2]
%   At         - Throat area [m^2]
%   Rcc        - Combustion chamber radius [m]
%   thickness  - Wall thickness [m]
%   Lcc        - just to plot CC
%   plots     = 1 if you want the plot outputs, = 0 if not
%
% Outputs:
%   y_patched  - 
%   r_nodes    - [x_nodes, y_nodes] 2xN
%   tangents   - 2xN : tangent vectors - tangents of the UPPER profile of
%                      the Rao Bell nozzle
%   normals    - 2xN : normal vectors - normals of the UPPER profile of the
%                      Rao bell Nozzle, UPWARD !
%   Lc_Rao     - Convergent length (Rao) [m]
%   Ld_Rao     - Divergent length (Rao) [m]
%   y_RAO_vec  - Radial coordinate vector of the nozzle profile [m]
%
% Notes:
%   - The Rao contour is approximated using a classical bell approach
%   - The method ensures similar divergence losses compared to the
%     equivalent conical nozzle.
%   - The geometry is continuous in slope at patching points.
%   - Includes 2D profile visualization and 3D surface reconstruction.
%
% References:
%   Rao, G.V.R. (1958), "Exhaust Nozzle Contour for Optimum Thrust"
%   F.Maggi, Slides - Non Reacting Flows
%  Fonti slide esercitazione:
%  G. P. Sutton and O. Biblarz. 2010. Rocket Propulsion Elements. John Wiley & Sons, Inc. Hoboken, NJ, USA.
%  C. Rogers. 2004. The Solid Rocket Motor – Part 4, Departures from Ideal Performance for Conical Nozzles and Bell Nozzles, Straight-Cut Throats and Rounded Throats. Tech-Series.

%
%--------------------------------------------------------------------------
% Author: Giorgio Vedovi, Sara Pino

% Input

% Input: alpha conico alpha [deg]
% Calcola Conic losses e conic sizes
dt = 2*sqrt(At/pi);         
de = 2*sqrt(Ae/pi);
Rt = dt/2;
Re = de/2;
Lc_conic = (Rcc - dt/2)/(tan(beta));
Ld_conic  = (Re - Rt)/(tan(alpha));
Ld_conic15 = (Re - Rt)/(tan(pi/12)); % 15° case
Lnozzle = Lc_conic + Ld_conic;

lambda = (1+cos(alpha))/2;

% calcola Rao che ha le stesse divergence losses (serve grafico parametrizzato)
% hp 
% N.B. since in the bell you impose that you want same beta and Rcc, you will 
% find a different convergent length, and won't find Lc_Rao = Lc_conic !
% via grafico
 % eta = %L wrt conical with 45°: find %L = eta s.t. same losses lambda
 % Having fixed eps, we iterate on eta = %L and find the iso-losses Rao
 % geometry:
 eta_isoloss = 0;
 differ = inf;
 % lambda_RAO_iso = 0;
 eps = Ae/At;
 for eta = 0.6 : 0.001 : 1
    alpha1 = atan((Re-Rt)/(eta*Ld_conic15));
    [theta_i, theta_e] = RaoCurves(eps, eta);
    th_i = deg2rad(theta_i);
    th_e = deg2rad(theta_e);
    lambda_RAO = 1/2 * (1 + cos((alpha1 + th_e)/2));
      if abs(lambda_RAO - lambda) < differ
          eta_isoloss = eta;
          % lambda_RAO_iso = lambda_RAO ;
          differ = abs(lambda_RAO - lambda);
          final_theta_i = th_i;
      end
 end


Ld_Rao = eta_isoloss*Ld_conic15;



% ----- BELL GEOMETRY ----------


% Rao approximated Bell Geometry
R1 = 1.5*Rt;
R2 = 0.382*Rt;

% Convergent-part circonference
y_C1 = @(x) R1 - sqrt (R1^2 - x.^2);
% Divergent-part circonference
y_C2 = @(x) R2 - sqrt (R2^2 - x.^2);


% Point 3 [Circ. 2 and parabola intersection]
x3 = R2*sin(th_i);
y3 = y_C2(x3);

% Point 4:
x4 = Ld_Rao;
y4 = Re - Rt;

% Dati di input (Esempio)
% x3, y3, x4, y4, th_i (in radianti)
m3 = tan(pi/2 - th_i); % La pendenza m3 è la tangente dell'angolo iniziale

% Costruzione della matrice dei coefficienti M
% riga 1: eq. nel punto 3
% riga 2: eq. nel punto 4
% riga 3: derivata nel punto 3
M = [y3^2,   y3,   1; ...
     y4^2,   y4,   1; ...
     2*y3,    1,   0];

% Vettore dei termini noti
Y_vec = [x3; x4; m3];

% Risoluzione del sistema lineare p = [a; b; c]
p = M \ Y_vec;

a = p(1);
b = p(2);
c = p(3);

parabola_params = [a;b;c];

x_par = @(x) a*x.^2 + b*x + c;
y_P = @(x) (-b+sqrt(b^2-4*a*(c-x)))/(2*a);
%plot(linspace(y3, y4, 1000), x_par(linspace(y3, y4, 1000)))
% Point 2
x2 = 0 ;
y2 = 0;

% % Point 1
% % tanBeta = (Rcc-y1)/(Lc_conic + x1)
x1 = -R1*sin(beta);
y1 = R1 - sqrt(R1^2 - x1^2);

% Line
m =  - tan(beta); 
q = y1 - m*x1;  

y_Line = @(x) m.*x + q;

% Point 0:
x0 = x1 - (Rcc - Rt - y1)*tan(beta);


% Patched nozzle profile function
y_patched = @(x) y_Line(x) .* (x >= x0 & x < x1) + ...
                 y_C1(x)  .* (x >= x1 & x < x2) + ...
                 y_C2(x)  .* (x >= x2 & x < x3) + ...
                 y_P(x)  .* (x >= x3 & x <= x4) + Rt;


x_vec = x0:0.0001:x4;
y_vec = y_patched(x_vec);

if plots == 1
figure; axis equal;
plot(x_vec,y_vec); hold on;
plot(x_vec,-y_vec); hold on;
end

Lc_Rao = (Rcc - Rt - y1)/tan(beta) - x1;

%% NODES GENERATION


l_Line = @(x1,x2) sqrt((x1-x2)^2 + (y_Line(x1)-y_Line(x2))^2);
l_C1   = @(x1,x2) R1*(atan((R1-y_C1(x2))/abs(x2)) - atan((R1-y_C1(x1))/abs(x1)));
l_C2   = @(x1,x2) R2* (asin(x2/R2) - asin(x1/R2));
l_P    = @(x1,x2) 1/(4*a).*( (asinh(2*a*y_P(x2)+b)+(2*a*y_P(x2)+b)*sqrt(1+(2*a*y_P(x2)+b)^2)) - (asinh(2*a*y_P(x1)+b)+(2*a*y_P(x1)+b)*sqrt(1+(2*a*y_P(x1)+b)^2)));

L1 = l_Line(x0,x1);
LC1 = l_C1(x1,x2);
LC2 = l_C2 (x2,x3);
L_P = l_P(x3,x4);

L =  L1 + LC1 + LC2 + L_P; % total length


dl = L/(N-1); % N nodes -> N-1 elements/curvilinear segments

% Line
N1 = 1e7; % discretize the l(x) very thin s.t. x(l) almost continuous
x_nodes1 = zeros(1,N);
x01_vec = linspace(x0, x1, N1);
l_Line_x = zeros(1,N1);

for j = 1 : N1
l_Line_x(j) = l_Line(x0,x01_vec(j));
end

j = 0;
s = 0;
while s + dl <= L1
     j = j + 1;
     if j ~=1
     s = s + dl;
     end
     [~,idx] = min(abs(l_Line_x - s));
     x_nodes1(j) = x01_vec(idx);
end
x_nodes1 = nonzeros(x_nodes1)'; % row vector

% the first curvilinear step in the next line will be shorter, in
% particular: % last s has overflown the arc
dl0 = dl - (L1 - s);

% Circonference 1

N2 = 1e7;
x_nodes2 = zeros(1,N2);
x12_vec = linspace(x1, x2, N2);
l_circ1_x = zeros(1,N2);

for j = 1 : N2
l_circ1_x(j) = l_C1(x1,x12_vec(j));
end

j = 0;
s = dl0;
while s + dl <= LC1
     j = j + 1;
     if j ~=1
     s = s + dl;
     end
     [~,idx] = min(abs(l_circ1_x - s));
     x_nodes2(j) = x12_vec(idx);
end
x_nodes2 = nonzeros(x_nodes2)';

% the first curvilinear step in the next line will be shorter, in
% particular:  last s has overflown the arc
dl0 = dl - (LC1-s); 

% ------- Circonference 2

N3 = 1e7;
x_nodes3 = zeros(1,N3);
x23_vec = linspace(x2, x3, N3);
l_circ2_x = zeros(1,N3);

for j = 1 : N3
l_circ2_x(j) = l_C2(x2,x23_vec(j));
end

j = 0;
s = dl0;
while s <= LC2 + 1e-12
    j = j + 1;
    [~,idx] = min(abs(l_circ2_x - s));
    x_nodes3(j) = x23_vec(idx);
    s = s + dl;
end
x_nodes3 = x_nodes3(1:j);
dl0 = s - LC2;

% ------- Parabola

N4 = 1e7;
x_nodes4 = zeros(1,N4);
x34_vec = linspace(x3, x4, N4);
l_par_x = zeros(1,N4);

for j = 1 : N4
l_par_x(j) = l_P(x3,x34_vec(j));
end

j = 0;
s = dl0;
while s + dl <= L_P
     j = j + 1;
     if j ~=1
     s = s + dl;
     end
     [~,idx] = min(abs(l_par_x - s));
     x_nodes4(j) = x34_vec(idx);
end
x_nodes4 = nonzeros(x_nodes4)';
if x_nodes4 ~= x4
    x_nodes4 = [x_nodes4, x4];
end
%x_nodes4(end)
%residual = (L_P - s); % as expected, the residual between the last s and the end of the nozzle is almost 0, i.e. I have arrived to the node = end of nozzle
%x_nodes4(end) - x4 % = 0. Ok.
% OBS.: the last curvilinear element will be slightly
% different due to non infinite discretization with which I advanced along
% the curves, accumulating errors at each passage from one curve to
% another).
% --------------

%figure;
% plot(xvec,y_patched(xvec));
x_nodes = [x_nodes1(1:end),x_nodes2(1:end),x_nodes3(1:end),x_nodes4(1:end)];
%size(x_nodes) % ok: 1xN

% PLOT NODES

y_nodes = y_patched(x_nodes);

if plots == 1
figure; axis equal;
plot(x_nodes,  y_nodes, '-ob'); hold on;
plot(x_nodes, -y_nodes, '-ob');
end

r_nodes = [x_nodes; y_nodes];


%% Analytical derivatives of the profile, for exact/analytical node extrusion direction
% spacing for discretization into adiabatic finite nozzle wall elements

dy_Line = @(x) m * ones(size(x));
dy_C1   = @(x) x ./ sqrt(R1^2 - x.^2);
dy_C2   = @(x) x ./ sqrt(R2^2 - x.^2);
dy_P    = @(x) 1 ./ sqrt(b^2 - 4*a*(c - x));  % derivata dell'inversa

dy_patched = @(x) dy_Line(x) .* (x >= x0 & x < x1) + ...
                  dy_C1(x)   .* (x >= x1 & x < x2) + ...
                  dy_C2(x)   .* (x >= x2 & x < x3) + ...
                  dy_P(x)    .* (x >= x3 & x <= x4);

% check: must be at least C1 by imposed bitangence at geometry switches:
%figure;
%plot(x_vec, dy_patched(x_vec)); % continuous, ok

% Tangent and normal unit vectors:
% Derivatives @ nodi x_nodes:
dy = dy_patched(x_nodes);          % 1xN
dx = ones(size(x_nodes));          % 1xN

% Norm-2
norms = sqrt(dx.^2 + dy.^2);      % 1xN

% Versori tangenti (2xN)
tangents = [dx ./ norms; dy ./ norms];

% Versori normali (2xN) — outward
normals = [-dy ./ norms; dx ./ norms];

%% 3D PLOT 
if plots == 1

theta = linspace(0, 2*pi, 200);
[X, TH] = meshgrid(x_vec, theta);


figure; hold on; grid on;

% =========================
% SUPERFICIE INTERNA
% =========================
R_int = abs(y_vec);

Y_int = R_int .* cos(TH);
Z_int = R_int .* sin(TH);

surf(X, Y_int, Z_int);

% =========================
% SUPERFICIE ESTERNA (OFFSET)
% =========================
R_ext = R_int + thickness;

Y_ext = R_ext .* cos(TH);
Z_ext = R_ext .* sin(TH);

surf(X, Y_ext, Z_ext);

% =========================
% CORONA INIZIALE (x = x0)
% =========================
theta_ring = linspace(0,2*pi,200);

Rin_0  = R_int(1);
Rout_0 = R_ext(1);

[TH0, RR0] = meshgrid(theta_ring, [Rin_0 Rout_0]);

X0 = x_vec(1) * ones(size(RR0));
Y0 = RR0 .* cos(TH0);
Z0 = RR0 .* sin(TH0);

surf(X0, Y0, Z0);

% =========================
% CORONA FINALE (x = x4)
% =========================
Rin_4  = R_int(end);
Rout_4 = R_ext(end);

[TH4, RR4] = meshgrid(theta_ring, [Rin_4 Rout_4]);

X4 = x_vec(end) * ones(size(RR4));
Y4 = RR4 .* cos(TH4);
Z4 = RR4 .* sin(TH4);

surf(X4, Y4, Z4);


% Combustion Chamber con spessore
x_cc = linspace(x_vec(1) - Lcc, x_vec(1), 100);
[XC, THC] = meshgrid(x_cc, theta);

% Superficie INTERNA (raggio Rcc)
R_in = Rcc * ones(size(XC));
surf(XC, R_in.*cos(THC), R_in.*sin(THC));

% Superficie ESTERNA (raggio Rcc + thickness)
R_out = (Rcc + thickness) * ones(size(XC));
surf(XC, R_out.*cos(THC), R_out.*sin(THC));

% Faccia anulare SINISTRA (fondo)
r_ring = linspace(Rcc, Rcc+thickness, 20);
[THb, RRb] = meshgrid(theta, r_ring);
Xb_left = (x_vec(1) - Lcc) * ones(size(RRb));
surf(Xb_left, RRb.*cos(THb), RRb.*sin(THb));

% Faccia anulare DESTRA (giunzione con nozzle)
Xb_right = x_vec(1) * ones(size(RRb));
surf(Xb_right, RRb.*cos(THb), RRb.*sin(THb));


% ---------------

% RENDERING

shading interp
axis equal
xlabel('x'); ylabel('y'); zlabel('z');

colormap(gray)

lighting phong
camlight headlight
camlight left
camlight right

material([0.3 0.6 1.0 50 1])




axis on
view(3)



full_length = x4 - x0;
end

%% Nested functions (Useful just for Rao approximation purpose)

function [theta_i, theta_e] = RaoCurves(eps, perc)

% Author: Sara Pino
% perc: value between 0 and 1, e.g. 0.5 for 50% length of L conv with alpha = 45°
% eps = Ae/At

% theta_ i , theta_e in [deg]

% data derived form graph [REFERENCE: SLIDE 74 - NON REACTING FLOWS]

x_60_e = [4.976; 10.685; 19.902; 29.527; 40.457; 49.674];
y_60_e = [19.426; 16.639; 14.344; 13.115; 12.541; 12.377];


x_70_e = [5.057; 9.951; 19.657; 29.608; 39.560; 49.918];
y_70_e = [16.013; 13.725; 11.520; 10.376; 9.804; 9.559];

x_80_e = [4.894; 9.788; 19.739; 29.853; 39.478; 49.837];
y_80_e = [13.154; 10.539; 8.742; 8.333; 7.843; 7.435];

x_90_e = [4.976; 10.114; 19.494; 29.853; 39.886; 49.592];
y_90_e = [10.866; 8.170; 7.026; 6.618; 6.291; 5.964];

x_100_e = [4.812; 9.869; 19.657; 29.853; 39.723; 49.592];
y_100_e = [8.252; 6.127; 5.229; 4.820; 4.575; 4.412];




x_100_i = [5.302; 9.788; 20.392; 30.261; 40.620; 49.592];
y_100_i = [20.000; 22.049; 25.492; 27.049; 28.361; 28.934];

x_90_i = [5.139; 10.277; 20.555; 31.321; 40.375; 49.755];
y_90_i = [21.066; 23.689; 27.295; 28.934; 29.672; 30.164];

x_80_i = [5.465; 10.033; 20.228; 30.914; 40.212; 49.511];
y_80_i = [22.951; 25.328; 29.098; 30.574; 31.311; 31.885];

x_70_i = [5.302; 10.196; 20.147; 30.098; 39.886; 49.592];
y_70_i = [25.738; 28.279; 31.803; 33.361; 34.180; 34.508];

x_60_i = [5.139; 10.767; 20.473; 30.914; 40.865; 49.347];
y_60_i = [28.443; 32.705; 35.246; 36.803; 37.623; 38.197];



th_i_60 = spline(x_60_i, y_60_i, eps);
th_i_70 = spline(x_70_i, y_70_i, eps);
th_i_80 = spline(x_80_i, y_80_i, eps);
th_i_90 = spline(x_90_i, y_90_i, eps);
th_i_100 = spline(x_100_i, y_100_i, eps);

th_e_60 = spline(x_60_e, y_60_e, eps);
th_e_70 = spline(x_70_e, y_70_e, eps);
th_e_80 = spline(x_80_e, y_80_e, eps);
th_e_90 = spline(x_90_e, y_90_e, eps);
th_e_100 = spline(x_100_e, y_100_e, eps);  


perc_vec = [0.60; 0.70; 0.80; 0.90; 1];
vect_th_i = [th_i_60; th_i_70; th_i_80; th_i_90; th_i_100];
vect_th_e = [th_e_60; th_e_70; th_e_80; th_e_90; th_e_100];


theta_i = spline(perc_vec, vect_th_i, perc);
theta_e = spline(perc_vec, vect_th_e, perc);

end



end

