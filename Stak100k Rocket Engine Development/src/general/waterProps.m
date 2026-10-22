function props = waterProps(T_C, P_bar)
 % CAREFUL INPUT IS IN CELSIUS AND BAR
T_ref = [10, 20, 30, 40, 50, 60, 70, 80];  
cp_ref = [4192, 4182, 4179, 4179, 4181, 4185, 4190, 4197]; 
mu_ref = [1307, 1002, 798, 653, 547, 467, 404, 354] * 1e-6; 
k_ref  = [0.580, 0.598, 0.615, 0.631, 0.644, 0.654, 0.663, 0.670]; 
rho_ref= [999.7, 998.2, 995.7, 992.2, 988.1, 983.2, 977.8, 971.8]; 
props.cp  = interp1(T_ref, cp_ref,  T_C, 'pchip', 'extrap');
props.mu  = interp1(T_ref, mu_ref,  T_C, 'pchip', 'extrap');
props.k   = interp1(T_ref, k_ref,   T_C, 'pchip', 'extrap');
props.rho = interp1(T_ref, rho_ref, T_C, 'pchip', 'extrap');
props.Pr  = props.mu * props.cp / props.k;
A = 5.20389; B = 1733.926; C = -39.485;
props.T_sat = (B / (A - log10(P_bar))) - C;% OUTPUT IN KELVIN
end