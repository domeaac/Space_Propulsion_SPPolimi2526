function K = BondCoat(T)

A = -6.3*1e-11; % W/(m K^5)
B = 2*1e-7;   % W/(m K^4)
C = -2*1e-4;  % W/(m K^3)
D = 1*1e-1;   % W/(m K^2)
E = -6.69 ;     % W/(m K^1)

K_BC = @(T) A.*T.^4 + B*T.^3 + C*T.^2 + D*T + E; % W/(mK)   % from "Thermophysical properties of Cold and Vacuum Plasma Spried CuCr..." model validity: 296 to 1223 K

K = K_BC(T);

end