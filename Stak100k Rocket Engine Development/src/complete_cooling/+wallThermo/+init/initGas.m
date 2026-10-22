function Gas = initGas(options)
% INITGAS Initializes the thermophysical properties of the exhaust gas.
%
% Inputs (All Optional / Name-Value Pairs):
%   options.cstar   - Characteristic velocity (m/s)
%   options.Prandtl - Prandtl number
%   options.P0      - Chamber stagnation pressure (Pa)
%   options.gamma   - Specific heat ratio
%   options.mu0     - Dynamic viscosity at stagnation temp (kg/(m*s))
%   options.cp      - Specific heat at constant pressure (J/(kg*K))
%
% Outputs:
%   Gas             - Struct containing all initialized gas properties.

    arguments
        options.cstar = 1426;
        options.Prandtl = 0.6;
        options.P0 = 70*1e5;
        options.gamma = 1.238;
        options.mu0 = 0.781*1e-4;
        options.cp = 1.9785;
    end 

    Gas = options;
    
end