%% CHAMBER DEFLATION SIMULATION
[t2,sol2, mass_rate2, Thrust2, subroutine_time_end_deflation] = deflation_RK4(grain, p_cut(end), timestep, tend/10, 5e5); 

warning('off','all')
tadditional = t2(2:subroutine_time_end_deflation/timestep);
soladditional = sol2(2:subroutine_time_end_deflation/timestep);
thrustadditional = Thrust2(2:subroutine_time_end_deflation/timestep);
warning('on','all')


t_cut = [t_cut, tadditional + t_cut(end)];
p_cut = [p_cut, soladditional];
Thrust_cut = [Thrust_cut; thrustadditional];

