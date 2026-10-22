function h_hotGas = convHeatCoefBartz(mu0, cp,Rc_t,cstar, Pr, P0,At,localArea, Twg, Tc, gamma, M)

    Dt = sqrt(4/pi*At); % throat diameter in meters 
    
    sigma = (0.5 .* (Twg ./ Tc) .* (1 + (gamma - 1) / 2 .* M.^2) + 0.5).^-0.68 .* (1 + (gamma - 1) / 2 .* M.^2).^-0.12;

    h = 0.026 / Dt^(0.2) * (Dt / Rc_t)^(0.1) * mu0^0.2 * cp / (Pr^0.6) * (P0/cstar)^(0.8) * (At / localArea)^0.9 .* sigma; 

    h_hotGas = h; 
    

end