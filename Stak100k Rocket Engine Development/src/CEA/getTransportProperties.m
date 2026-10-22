function [ceaout] = getTransportProperties(grain)

ceaout = wrp_rkt_AP_HTPB_tran(grain.design.p0/1e5, grain.design.p0_over_pe, grain.temperature, ...
    grain.formulation.AP_mass_fraction, grain.formulation.HTPB_mass_fraction);

end