function T_func = getParabolicProfile(x_v,Tmax, a)
% GETPARABOLICPROFILE Generates a parabolic function handle based on vertex conditions.
%
% Inputs:
%   x1, x2 - Coordinates used to define the vertex location (x_v = x2 - x1)
%   Tmax   - The maximum value of the parabola
%   a      - (Optional) Concavity/spread factor. Must be negative. Default is -1.
%
% Output:
%   T_func - Anonymous function handle that evaluates the parabola

    arguments
        
        x_v (1,1) double
        Tmax (1,1) double
        a (1,1) double  = 1
    end
    
   
    
    % Create and return the anonymous function in vertex form
    T_func = @(x) a * (x - x_v).^2 + Tmax;

end