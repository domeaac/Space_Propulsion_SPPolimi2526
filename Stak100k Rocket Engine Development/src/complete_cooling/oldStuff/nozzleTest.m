alpha = deg2rad(10); % divergent angle (conic nozzle)
beta = deg2rad(45); % convergent angle (conic nozzle)
At = 0.0089; % m^2
eps = 8;
Ae = 8*At;
Rcc = 0.432;
thickness = 0.5;

% Input

% Input: alpha conico alpha [deg]
% Calcola Conic losses e conic sizes
dt = 4*sqrt(At/pi);         
de = 4*sqrt(Ae/pi);
Rt = dt/2;
Re = de/2;
Lc_conic = (Rcc - dt/2)/(tan(beta));
Ld_conic  = (Re - Rt)/(tan(alpha));
Ld_conic45 = (Re - Rt)/(tan(pi/4)); % 45° case
Lnozzle = Lc_conic + Ld_conic;

lambda = (1+cos(alpha))/2;

% calcola Rao che ha le stesse divergence losses (serve grafico parametrizzato)
% hp 
Lc_Rao = Lc_conic; % supposed equal in Rao
% via grafico
eta = 0.7; % %L wrt conical with 45° (!)
th_e = deg2rad(15);
th_i = deg2rad(27.5);
% solve eqz. s.t. same losses.
Ld_RAO = eta*Ld_conic45;

% ----- BELL GEOMETRY ----------


% Rao approximated Bell Geometry
R1 = 1.5*Rt;
R2 = 0.382*Rt;

% Convergent-part circonference
y_C1 = @(x) R1 - sqrt (R1^2 - x.^2);
% Divergent-part circonference
y_C2 = @(x) R2 - sqrt (R2^2 - x.^2);


% Point 3 [Circ. 2 and parabola intersection]
x3 =R2*sin(th_i);
y3 = y_C2(x3);

% Parabola
m = tan(pi/2-th_i);
n = tan(pi/2-th_e);

a = (n-m)/(2*(Re-Rt)-2*y3);
b = m - 2*a*y3;
c = x3-a*y3^2 + (2*a*y3-m)*y3;

y_P = @(x) (-b+sqrt(b^2-4*a*(c-x)))/(2*a);

% Point 2
x2 = 0 ;
y2 = 0;

% % Point 1
% % tanBeta = (Rcc-y1)/(Lc_conic + x1)

x1 = -R1*sin(beta);
y1 = R1 - sqrt(R1^2 - x1^2);



% Line
m = -tan(beta);
q = y1 - m*x1;
y_Line = @(x) m.*x + q;

% Point 4:
x4 = Ld_RAO;
y4 = Re - Rt;



% Point 0:
x0 = -Lc_Rao;


% Patched nozzle profile function
y_patched = @(x) y_Line(x) .* (x >= x0 & x < x1) + ...
                 y_C1(x)  .* (x >= x1 & x < x2) + ...
                 y_C2(x)  .* (x >= x2 & x < x3) + ...
                 y_P(x)  .* (x >= x3 & x <= x4) + Rt;


x_vec = x0:0.0001:x4;
y_vec = y_patched(x_vec);

figure("WindowStyle", "docked"); axis equal;

plot(x_vec,y_vec); hold on;
plot(x_vec,-y_vec); hold on;



theta = linspace(0, 2*pi, 200);  

[X, TH] = meshgrid(x_vec, theta);

R = abs(y_vec);   % positive radius check

Y = R .* cos(TH);
Z = R .* sin(TH);

figure;
surf(X, Y, Z);

% Rendering 
shading interp
axis equal
xlabel('x'); ylabel('y'); zlabel('z');

% Color gray ~ metal
colormap(gray)

% Shadows and illumination 
lighting phong
camlight headlight
camlight left
camlight right

% Materiale cromo
material([0.3 0.6 1.0 50 1])  
% [ambient diffuse specular spec_exp spec_color]

% Visualization settings
%set(gcf, 'Color', 'k')     % black background
%set(gca, 'Color', 'k')     % black axes background
axis on
%axis off

