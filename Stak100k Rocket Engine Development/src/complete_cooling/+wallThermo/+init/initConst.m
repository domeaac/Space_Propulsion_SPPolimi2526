function Const = initConst(options)
% INITCONST Initializes the constant parameters for the nozzle solver.
%
% Inputs:
%   options.TotalTemperature - Chamber stagnation temperature (K) [Default: 3000]
%   options.NozzleCurve      - Function handle for the nozzle profile y = g(x)
%   options.At               - Geometric throat area (m^2)
%   options.x_throat         - Axial coordinate of the throat (m)
%   options.Rc_t             - Radius of curvature at the throat (m)
%
% Outputs:
%   Const                    - Struct containing all initialized constants.

    arguments
        options.TotalTemperature (1,1) double = 3000;
        options.NozzleCurve = -1;
        options.x_throat = NaN;
        options.Rc_t = [];
    end 

  
    Const = options;
    Const.At = pi * options.NozzleCurve(options.x_throat)^2; 
    
end