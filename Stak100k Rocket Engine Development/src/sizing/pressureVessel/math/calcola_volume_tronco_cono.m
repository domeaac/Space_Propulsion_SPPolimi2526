function V = calcola_volume_tronco_cono(D1, theta_deg, L)
%CALCOLA_VOLUME_TRONCO_CONO Calcola il volume interno di un tronco di cono.
%
%   V = calcola_volume_tronco_cono(D1, theta_deg, L)
%
%   INPUT:
%   D1        - Diametro interno della base minore [m]
%   theta_deg - Angolo di inclinazione della parete rispetto all'orizzontale [gradi]
%   L         - Lunghezza assiale (altezza) del tronco di cono [m]
%
%   OUTPUT:
%   V         - Volume interno del tronco di cono [m^3]

R1 = D1 / 2;
R2 = R1 + L * cotd(theta_deg);
V = (pi * L / 3) * (R1^2 + R2^2 + R1 * R2);

end