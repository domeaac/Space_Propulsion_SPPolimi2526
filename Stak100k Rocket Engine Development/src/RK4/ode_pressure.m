function ddtpc=ode_pressure(t, p_c,  param , A_tr , rho_prop, rb, burning_area, inner_grain_volume)

%FUNCTION FOR INTEGRATION ddtP_c=f(p_c, t)

%-------------------------INPUT--------------------------------------

% t -> time [s]
% p_c -> chamber pressure [Pa]
% param -> output from cea at every timestep [struct from CEA]
% A_tr -> throat area for given grain fomrulation [m^2]
% rho_prop -> solid grain density [kg/m^3]
% rb (T , p_c) -> burning rate as function handle , p_c is in Mpa , rb is in m/sec

%-------------------------OUTPUT-------------------------------------

% ddtpc -> timederivative of pressure in Pa/s

%--------------------------------------------------------------------

scale=1e5; %1bar *scale= Pa
R = 8314; %gas constant 
rb = rb/1e3;

%gas density in CHAMBER 

rho_ch=param.output.froz.density(1); 

%c_star in EXIT 
c_star=param.output.froz.cstar(end); 

%massa molare in CHAMBER
M_w=param.output.froz.mw(1); 
% si utilizza massa con condensati, piu accurato di massa solo gas essendo
% l' equazione che segue un bilancio di massa.
% Il software NASA CEA ipotizza che la relazione gas perfetto sia valida anche
% in presenza di condensati (si faccia riferimento a NASA RP‑1311 inizio
% capitolo 2.2). La scelta di usare mw condensati, dunque, non aggiunge ipotesi aggiuntive al
% processo di design ed assicura la coerenza tra i vari ternini del
% bilancio di massa che segue. Scegliere mw senza condensati
% causerebbe una discordanza tra le approx usate nei due termini: cstar
% infatti è calcolato includendo la massa dei condensati (si faccia
% riferimento a NASA RP‑1311)
% Questa scelta, inoltre, permette di mantenere la concordanza tra il
% modello di integrazione equilibrium ed il modello di integrazione frozen;
% per quest' ultimo si è infatti ricevuta l' esplicita indicazione di utilizzare mw
% contenente condensati.
 
%temp in CHAMBER
temp=param.output.froz.temperature(1); 

ddtpc = ( (rho_prop-rho_ch).* burning_area.* rb - p_c.*A_tr./c_star) .* (inner_grain_volume/(R/M_w*temp)).^(-1);

% rho_prop
% rho_ch
% burning_area
% rb
% p_c/1e6
% A_tr
% c_star




%---------------------------REAL TIME OUTPUT-------------------------------
fprintf('pressure : %f Bar \n', p_c/scale)
fprintf('pressure Gradient : %f bar/s \n', ddtpc/scale) 
%fprintf('burning rate %f m/s \n', rb) 
% fprintf('BBB %f \n', ( (rho_prop-rho_ch).* burning_area.* rb - p_c.*A_tr./c_star)) 
% fprintf('AAA %f \n', (inner_grain_volume/(R/M_w*temp)).^(-1)) 
fprintf('time : %f s \n', t) 


end 