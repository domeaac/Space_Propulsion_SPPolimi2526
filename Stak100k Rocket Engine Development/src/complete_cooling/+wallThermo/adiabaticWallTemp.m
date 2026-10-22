function T_aw = adiabaticWallTemp(gamma, Pr,MachNumber, T_ref)
% ADIABATICWALLTEMP Computes the local adiabatic wall temperature of the nozzle.

    % Compute turbulent recovery factor (T_aw / T_0)
    R = (1 + Pr^(1/3) * (gamma-1)/2 * MachNumber^2) / (1 + (gamma-1)/2 * MachNumber^2); 

    T_aw = R * T_ref; 

end