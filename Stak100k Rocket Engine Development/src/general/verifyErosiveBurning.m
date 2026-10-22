function [r, r0, re, ratio] = verifyErosiveBurning(grain)

    r0 = grain.r_b(grain.design.p0/1e5);
    ceaout = wrp_rkt_AP_HTPB_tran(grain.design.p0/1e5, grain.design.p0_over_pe, grain.temperature, ...
    grain.formulation.AP_mass_fraction, grain.formulation.HTPB_mass_fraction);
    
    pr = ceaout.output.froz.prandtl.froz(1);
    mu = ceaout.output.froz.viscosity(1) * 1e-4; %1e-6 ???
    cp = ceaout.output.froz.cp(1) * 1e3;   
    rhob = grain.formulation.rho_av(1);
    T1 = ceaout.output.froz.temperature(1);
    Tb = grain.temperature(1);
    beta = 53;
    D = (4*pi*grain.Ri0^2)/(2*pi*grain.Ri0);
    G = ((grain.Re0^2 - grain.Ri0^2) + grain.L0*grain.Ri0*2*pi) * (grain.r_b(grain.design.p0/1e5)/1e3)*grain.formulation.rho_av ;

    cs = 1200; % 
    Ts = 800; % from Kubota [K]
    T2 = T1;  % 

    alpha = (0.0288*cp*(mu^0.2)*(pr^(-2/3))/(rhob*cs))*(T1-Ts)/(T2-Tb);
    re = alpha * (G^0.8) * (D^(-0.2)) * exp(-beta*r0/G);


    re = re*1e3;
    ratio = re/r0;
    r = r0 + re;

end