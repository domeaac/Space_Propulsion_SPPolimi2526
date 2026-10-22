function [M, prandtl, gamma, cp, mu, cstar, temperature] = getCEAdata4bera(suborsup, aeat_span, grain)

M = [];
prandtl = [];
gamma = [];
cp = [];
mu = [];
cstar = [];
temperature = [];

maxsiz = 10;

for i = 0:ceil(max(size(aeat_span))/maxsiz)-1
    if i == ceil(max(size(aeat_span))/maxsiz)-1
      aeat_vec = aeat_span((i*maxsiz + 1):end) %#ok<*FNCOLND>
    else
      aeat_vec = aeat_span((i*maxsiz + 1):(i*maxsiz + maxsiz))
    end

    [cea_output] = wrp_rkt_AP_HTPB_froz_tran_stations(grain.design.p0/1e5, aeat_vec, suborsup, grain.temperature, grain.formulation.AP_mass_fraction, grain.formulation.HTPB_mass_fraction);
   
    %cea_output.output.froz.pressure(3:end)
    M = [M; cea_output.output.froz.mach(3:end)];
    prandtl = [prandtl; cea_output.output.froz.prandtl.froz(3:end)];
    gamma = [gamma; cea_output.output.froz.gamma(3:end)];
    cp = [cp; cea_output.output.froz.cp(3:end)];
    mu = [mu; cea_output.output.froz.viscosity(3:end)];
    cstar = [cstar; cea_output.output.froz.cstar(3:end)];
    temperature = [temperature; cea_output.output.froz.temperature(3:end)];

end

end

