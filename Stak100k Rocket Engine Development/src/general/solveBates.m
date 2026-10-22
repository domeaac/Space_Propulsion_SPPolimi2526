function [BatesMotorStruct] = solveBates(grain, rb, massFlowRate,burningTime,rho_prop,opt)
    arguments
        grain
        rb
        massFlowRate
        burningTime
        rho_prop
        opt.x0  (1,3) double = [0.4 0.1 0.3]; 
    end  
    
    ceaout = grain.initial_CEA_run.output.froz;

    % assemble system of equations 
    totalMass = massFlowRate*burningTime; 
    
    Ab = @(h0 , d0 , dext , t) 2*pi/4 * (dext^2 - (d0 + 2*rb*t)^2 ) + pi*(d0+2*rb*t)*(h0 - 2*rb*t); 
    Vol = @(h0 , d0 , dext , t) pi/4 * (dext^2 - (d0 + 2*rb*t)^2 )*(h0 - 2*rb*t); 
     
    % eq 1 initial Area = final Area 
    f1 = @(h0, d0, dext) Ab(h0, d0, dext, 0)-Ab(h0,d0,dext,burningTime); 

    % eq 2 initial massflowrate 
    % f2 = @(h0, d0, dext) Ab(h0, d0, dext, 0)*rho_prop*rb - massFlowRate; 
    % eq 2.2 initial mass flowrate from unsteady equation mass balance
    %rho_ch= ceaout.density(1); 
    f2 = @(h0, d0, dext) (grain.formulation.rho_av-ceaout.density(1))*Ab(h0,d0,dext,0)*(grain.design.rb/1e3) - grain.design.p0*grain.at/grain.cstar;

    % eq 3 totalMass 
    f3 = @(h0, d0, dext) Vol(h0, d0, dext, 0)*rho_prop - totalMass; 


    % MAGGI SAID
    % f1 = Ab0 == Abf
    % f2 = P(t = 0) == Pdesign == (a c* (Ab/At) rho_prop)^1/(1-n)  
    % f3 = T(t = 0) == Tdesign == mdot ve == mdot Isp ve 
    % P = 70 bar is MEOP 
    % enforcing global quantities related equation lead to an error

    f = @(x) [f1(x(1), x(2), x(3)); f2(x(1), x(2), x(3)); f3(x(1), x(2), x(3))]; % x = [h0 d0 dext]

    % solve eq 
    x0 = opt.x0;
    options = optimoptions('fsolve', ...
    'Display', 'none', ...                 % Change to 'iter' if you need to debug convergence
    'Algorithm', 'trust-region-dogleg', ...% The default and usually best root-finder
    'FunctionTolerance', 1e-12, ...        % How close the equations (f1, f2, f3) must be to exactly 0
    'StepTolerance', 1e-12, ...            % Minimum allowed change in your variables (h, D_int, D_ext)
    'MaxFunctionEvaluations', 5000, ...    % Prevents the solver from giving up too early
    'MaxIterations', 2000);                % Maximum number of algorithmic steps
    
    sol = fsolve(f, x0,options); 

    BatesMotorStruct.h = sol(1); 
    BatesMotorStruct.D_int = sol(2); 
    BatesMotorStruct.D_ext = sol(3);
    
end