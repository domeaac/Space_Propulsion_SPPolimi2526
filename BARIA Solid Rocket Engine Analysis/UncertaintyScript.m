% sigma_a =
% sigma_n =
% Assuming all variables for all batches to be indipendent and identically
% distributed (same history of production), we can trace a regression line
% using all 27 points together:

X = log(p_val);                                                        % conversione valori di pressione in scala logaritmica
N = length(p_val);                                                     % lunghezza vettore delle pressioni
Y = log(r_val);                                                        % conversione valori di rb in scala logaritmica

delta = N * sum(X.^2)-(sum(X))^2;                                   % denominatore per il calcolo di q e m
m = ( N*sum(X.*Y) - sum(X)*sum(Y) ) / delta;
% Mean value of pressure exponent:
n = m;                                                              % dal passaggio sopra
q = ( sum(X.^2)*sum(Y) - sum(X)*sum(X.*Y)) / delta;                 % intercetta della retta Y = mX + q con m = n

Y_eval = m .* X + q;                                                % valori valutati tramite equazione della retta appena calcolata
sig = ( ( sum( ( Y - Y_eval ).^2) ) / (N-2) )^0.5;                  % scarto quadratico medio dei valori reali da quelli della retta

Inc_m = sig * (N/delta)^0.5;                                        % incertezza del coeff "m" (che coincide con quella di "n")
% Estimation of Standard deviation of pre-exponential coefficient
Inc_n = Inc_m;     % m=n ->  Inc_n = Inc_m                                                  % dal passaggio sopra
Inc_q = sig * ( sum( X.^2 ) / delta )^0.5;                          % incertezza del coeff "q"

% Estimation of Mean value of pre-exponential coefficient
a = exp(q);                                                         % calcolo del valore "a"
a_new_up = a * exp(Inc_q);                                          % calcolo della parte superiore l'incertezza di "a"
a_new_down = a / exp(Inc_q);                                        % calcolo della parte inferiore per l'incertezza di "a"
% Estimation of Standard deviation of pre-exponential coefficient
Inc_a = ( a_new_up - a_new_down) / 2;                               % calcolo di "a"


M_Y = mean(Y);                                                               % calcolo della media di Y
dev_reg = sum ( (Y_eval - M_Y).^2 );                                         % calcolo della devianza di regressione
dev_tot = sum ( (Y - M_Y).^2 );                                              % calcolo della devianza totale
R2 = dev_reg / dev_tot  % close to 1 -> data have strong linear dependance   % calcolo del coefficiente R^2

