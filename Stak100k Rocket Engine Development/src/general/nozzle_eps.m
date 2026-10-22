function eps = nozzle_eps(Pe_P0, gamma)
% NOZZLE_EPS  Calcola il rapporto di espansione epsilon = Ae/At
%
%   eps = nozzle_eps(Pe_P0, gamma)
%
%   Input:
%       Pe_P0  - rapporto Pe/P0 (deve essere 0 < Pe_P0 < 1)
%       gamma  - rapporto dei calori specifici (default: 1.4)
%
%   Output:
%       eps    - rapporto di espansione Ae/At

    if nargin < 2
        gamma = 1.4;
    end

    g = gamma;

    term1 = ((g + 1) / 2) .^ (1 / (g - 1));
    term2 = Pe_P0 .^ (1 / g);
    inner = 1 - Pe_P0 .^ ((g - 1) / g);
    term3 = sqrt( ((g + 1) / (g - 1)) .* inner );

    inv_eps = term1 .* term2 .* term3;
    eps = 1 ./ inv_eps;

end