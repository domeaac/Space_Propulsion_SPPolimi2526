function [t,sol, mass_rate, Thrust, time_end_combustion] = reduced_RK_2_frozen_Z(grain, x0, timestep, tend)

% OUTPUTS:
% t -> time vector
% sol -> p_c(t)
% mass_rate -> mass_rate(t)
% Thrust -> Thrust(t)


%-----------------------SELECT CURRENT FORMULATION PROPERTIES-------------
wrapper_CEA=str2func('wrp_rkt_AP_HTPB_imposing_aeat'); %wrapper selezionato


%initialise time vector
t=0:timestep:tend;
N = length(t);

%pralloco vettore soluzione
sol=zeros(1,N);
sol(1)=x0;

mass_rate = [];
Thrust = [];
aus_vect = [];

flag_combustione_esaurita = false;

time_end_combustion = 0;

%calcolo K per ogni step e integro
for i=1:N-1
    if flag_combustione_esaurita == false

        %PARAMETRI PRIMO STEP RK
        %calcolo rb
        rb_firstStep = grain.r_b(sol(i)/1e5);
        %RUN CEA FOR PARAMETERS USED IN INTEGRATION
        param_firstStep=wrapper_CEA(sol(i)/1e5, grain.aeat, grain.temperature, grain.formulation.AP_mass_fraction, grain.formulation.HTPB_mass_fraction); %runno cea con pressione, aeat, Tgrano
        %calcolo parametri geometrici PRIMO STEP
        grain_aux = grain;
        burning_area = BurningArea(grain_aux);
        inner_grain_volume = ccVolume(grain_aux);
        %PRIMO STEP RK
        k1=ode_pressure(t(i), sol(i),  param_firstStep , grain.at , grain.formulation.rho_av, rb_firstStep, burning_area, inner_grain_volume);


        fprintf('Residual length = %f.3  [m] \n', grain.x(4)-grain.x(1))
        fprintf('Residual height = %f.3  [m] \n', grain.y(2))

        %PARAMETRI SECONDO STEP RK
        %calcolo rb
        rb_secondStep = grain.r_b((sol(i) + timestep*k1/2)/1e5);
        %RUN CEA FOR PARAMETERS USED IN INTEGRATION
        param_secondStep=param_firstStep;
        %wrapper_CEA((sol(i) + timestep*k1/2)/1e5, grain.aeat, grain.temperature, grain.formulation.AP_mass_fraction, grain.formulation.HTPB_mass_fraction); %runno cea con pressione, aeat, Tgrano
        %calcolo parametri geometrici SECONDO STEP
        grain_aux = applyRegression(grain, timestep/2, rb_secondStep);
        burning_area = BurningArea(grain_aux);
        inner_grain_volume = ccVolume(grain_aux);
        %SECONDO STEP RK
        k2=ode_pressure(t(i) + timestep/2, sol(i) + timestep*k1/2,  param_secondStep, grain.at , grain.formulation.rho_av, rb_secondStep, burning_area, inner_grain_volume);


        fprintf('Residual length = %f.3  [m] \n', grain.x(4)-grain.x(1))
        fprintf('Residual height = %f.3  [m] \n', grain.y(2))

        %PARAMETRI TERZO STEP RK
        %calcolo rb
        rb_thirdStep = grain.r_b((sol(i) + timestep*k2/2)/1e5);
        %RUN CEA FOR PARAMETERS USED IN INTEGRATION
        param_thirdStep=param_secondStep;
        %wrapper_CEA((sol(i) + timestep*k2/2)/1e5, grain.aeat, grain.temperature, grain.formulation.AP_mass_fraction, grain.formulation.HTPB_mass_fraction); %runno cea con pressione, aeat, Tgrano
        %calcolo parametri geometrici TERZO STEP
        grain_aux = applyRegression(grain, timestep/2, rb_thirdStep);
        burning_area = BurningArea(grain_aux);
        inner_grain_volume = ccVolume(grain_aux);
        %TERZO STEP RK
        k3=ode_pressure(t(i) + timestep/2, sol(i) + timestep*k2/2,  param_thirdStep, grain.at , grain.formulation.rho_av, rb_thirdStep, burning_area, inner_grain_volume);


        fprintf('Residual length = %f.3  [m] \n', grain.x(4)-grain.x(1))
        fprintf('Residual height = %f.3  [m] \n', grain.y(2))

        %PARAMETRI QUARTO STEP RK
        %calcolo rb
        rb_fourthStep = grain.r_b((sol(i) + timestep*k3)/1e5);
        %RUN CEA FOR PARAMETERS USED IN INTEGRATION
        param_fourthStep=param_secondStep;
        %wrapper_CEA(( sol(i) + timestep*k3)/1e5, grain.aeat, grain.temperature, grain.formulation.AP_mass_fraction, grain.formulation.HTPB_mass_fraction); %runno cea con pressione, aeat, Tgrano
        %calcolo parametri geometrici QUARTO STEP
        grain_aux = applyRegression(grain, timestep, rb_fourthStep);
        burning_area = BurningArea(grain_aux);
        inner_grain_volume = ccVolume(grain_aux);
        %QUARTO STEP RK
        k4=ode_pressure(t(i) + timestep,   sol(i) + timestep*k3,    param_fourthStep, grain.at , grain.formulation.rho_av, rb_fourthStep, burning_area, inner_grain_volume);


        %AGGIORNO SOLUZIONE
        sol(i+1)=sol(i) + timestep/6 * (k1 + 2*k2 + 2*k3 + k4);

        %AGGIORNO GEOMETRIA CON BURNING RATE DI STEP 1 (calcolato alla pressione P(i))
        [grain, flag_combustione_esaurita] = applyRegression(grain,timestep,rb_firstStep);

        %AGGIORNO TRAZIONE E PORTATA MASSICA CON DATI DI STEP 1 (cioè calcolati alla pressione P(i))
        Thrust=[Thrust; grain.at*sol(i)*param_firstStep.output.froz.cf(3)*grain.LossesCoeff + (param_firstStep.output.froz.pressure(3)*1e5 - 101325)*grain.aeat*grain.at];
        aus_vect = [aus_vect; param_firstStep.output.froz.cf(3).*grain.at*sol(i)];
        mass_rate=[mass_rate;  param_firstStep.output.froz.pressure(1)*1e5*grain.at/param_firstStep.output.froz.cstar(end)];

        fprintf('Residual length = %f.3  [m] \n', grain.x(4)-grain.x(1))
        fprintf('Residual height = %f.3  [m] \n', grain.y(2))

        time_end_combustion = t(i);
    else
        sol(i) = -1;

    end
end

if flag_combustione_esaurita == false
%Salvo Trazione , portata massica per t=end
param = wrapper_CEA(sol(end)./1e5 , grain.aeat , grain.temperature, grain.formulation.AP_mass_fraction, grain.formulation.HTPB_mass_fraction);

p_c=param.output.froz.pressure(1)*1e5;
c_star=param.output.froz.cstar(end);

%Salvo Trazione , portata massica per ogni iterazione CEA
mass_rate=[mass_rate;  p_c*grain.at/c_star];
Thrust=[Thrust; grain.at*sol(i)*param.output.froz.cf(3)*grain.LossesCoeff + (param.output.froz.pressure(3)*1e5 - 101325)*grain.aeat*grain.at];
end


end
