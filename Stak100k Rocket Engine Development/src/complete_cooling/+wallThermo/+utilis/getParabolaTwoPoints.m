function P_func = getParabolaTwoPoints(x1, y1, x2, y2, a, options)
% GETPARABOLATWOPOINTS Generates a parabola passing through two points.
%
% Inputs:
%   x1, y1 - Coordinates of the first point
%   x2, y2 - Coordinates of the second point
%   a      - (Optional) The leading coefficient (concavity). Default is 1.
%   options.ZeroDerivativeAtX1 - (Optional logical) If true, forces a 
%                                horizontal tangent at (x1, y1). This 
%                                overrides the provided 'a' value.
%
% Output:
%   P_func - Anonymous function handle evaluating the parabola

    arguments
        x1 (1,1) double
        y1 (1,1) double
        x2 (1,1) double
        y2 (1,1) double
        a  (1,1) double = 1
        options.ZeroDerivativeAtX1 (1,1) logical = false
    end
    
    if x1 == x2
        error('getParabolaTwoPoints:InvalidInput', 'x1 and x2 cannot be the same value.');
    end

    if options.ZeroDerivativeAtX1
        % Condition 1: y'(x1) = 0  =>  2*a*x1 + b = 0  =>  b = -2*a*x1
        % Condition 2: Passes through (x1, y1)  =>  c = y1 + a*x1^2
        % Condition 3: Passes through (x2, y2)
        % Substituting b and c into Condition 3 yields the exact value for a:
        
        a = (y2 - y1) / (x2 - x1)^2;
        b = -2 * a * x1;
        c = y1 + a * x1^2;
        
    else
        % Solve the linear system for b and c using the prescribed 'a'
        M = [x1, 1; x2, 1];
        RHS = [y1 - a * x1^2; y2 - a * x2^2];
        
        coeffs = M \ RHS;
        b = coeffs(1);
        c = coeffs(2);
    end
    
    % Return the anonymous function
    P_func = @(x) a * x.^2 + b * x + c;
end