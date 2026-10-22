function Solution = solverColdSide(Mesh, k_layers, grain, options)
% SOLVERCOLDSIDE Computes steady-state 1D transverse thermal profiles for nozzle wall segments given a cold-side boundary condition.
%
% Inputs:
%   Mesh       - Struct containing geometric mesh data (curveNodes, nLayers, layerThickness).
%   k_layers   - 1D array of thermal conductivities for each wall layer [W/(m*K)].
%   sigma      - Bartz equation convective heat transfer correction factor [-].
%   grain      - Struct containing propellant and chamber properties.
%   options    - Struct containing optional function handles (BcTColdFuncHandle).
%
% Outputs:
%   Solution   - Struct containing the Mesh and evaluated thermal arrays.

    arguments 
        Mesh
        k_layers
        grain
        options.BcTColdFuncHandle = []; 
    end 

    Rc_t = Mesh.CurvatureRadiusThroat;
    At = Mesh.At;
    ThroatNodalPos = Mesh.throatPos;
    x_nodal_Throat = ThroatNodalPos(1);

    if Mesh.nLayers ~= length(k_layers)
        error("solverColdSide:IllegalLayer", "Number of layers must match the k_layers array length.");
    end

    x_points = Mesh.curveNodes(1,:);
    y_points = Mesh.curveNodes(2,:);

    x_midPanel = (x_points(1:end-1) + x_points(2:end)) / 2;
    y_midPanel = (y_points(1:end-1) + y_points(2:end)) / 2;
    nMidPoints = size(x_midPanel(:)', 2); 
    MidPointCoords = [x_midPanel(:)'; y_midPanel(:)']';

    if isempty(options.BcTColdFuncHandle)
        error("solverColdSide:error","No function handle prescribed for cold wall temperature")
    else 
        T_ColdWall = options.BcTColdFuncHandle(x_midPanel(:)); 
    end     

    if strcmp(Mesh.region, "sub-sup")
        mask_sub = x_midPanel <= x_nodal_Throat;
        x_pan_sub = x_midPanel(mask_sub);
        y_pan_sub = y_midPanel(mask_sub);
        Ax_pan_sub = pi * y_pan_sub.^2;
        AxAt_sub = Ax_pan_sub ./ At;

        mask_sup = x_midPanel > x_nodal_Throat;
        x_pan_sup = x_midPanel(mask_sup);
        y_pan_sup = y_midPanel(mask_sup);
        Ax_pan_sup = pi * y_pan_sup.^2 ;
        AxAt_sup = Ax_pan_sup ./ At;
        
        [M_sub, prandtl_sub, gamma_sub, cp_sub, mu_sub, cstar_sub, T_sub] = CEAThermo.CeaThermoSolver('sub', AxAt_sub, grain);
        [M_sup, prandtl_sup, gamma_sup, cp_sup, mu_sup, cstar_sup, T_sup] = CEAThermo.CeaThermoSolver('sup', AxAt_sup, grain);

        M = [M_sub; M_sup];
        prandtl = [prandtl_sub; prandtl_sup];
        gamma = [gamma_sub; gamma_sup];
        cp = [cp_sub; cp_sup];
        mu = [mu_sub; mu_sup];
        cstar = [cstar_sub; cstar_sup];
        T_cea = [T_sub; T_sup];
        AxAt = [AxAt_sub, AxAt_sup]';
        Ax = [Ax_pan_sub, Ax_pan_sup]';

    elseif strcmp(Mesh.region, "sub")
        x_pan_sub = x_midPanel; 
        y_pan_sub = y_midPanel;
        Ax_pan_sub = pi * y_pan_sub.^2 ;
        AxAt_sub = Ax_pan_sub ./ At;

        [M_sub, prandtl_sub, gamma_sub, cp_sub, mu_sub, cstar_sub, T_sub] = CEAThermo.CeaThermoSolver('sub', AxAt_sub, grain);

        M = M_sub;
        prandtl = prandtl_sub;
        gamma = gamma_sub;
        cp = cp_sub;
        mu = mu_sub;
        cstar = cstar_sub;
        T_cea = T_sub;
        AxAt = AxAt_sub';
        Ax = Ax_pan_sub';

    elseif strcmp(Mesh.region, "sup")
        x_pan_sup = x_midPanel;
        y_pan_sup = y_midPanel;
        Ax_pan_sup = pi * y_pan_sup.^2;
        AxAt_sup = Ax_pan_sup ./ At;

        [M_sup, prandtl_sup, gamma_sup, cp_sup, mu_sup, cstar_sup, T_sup] = CEAThermo.CeaThermoSolver('sup', AxAt_sup, grain);

        M = M_sup;
        prandtl = prandtl_sup;
        gamma =  gamma_sup;
        cp = cp_sup;
        mu = mu_sup;
        cstar = cstar_sup;
        T_cea = T_sup;
        AxAt = AxAt_sup';
        Ax = Ax_pan_sup';
    end

    Panels = [];
    N_MIDPOINTS = Mesh.nCurveNodes - 1;
    N_TEMPERATURES = Mesh.nLayers + 2;
    TemperatureProfile = zeros(N_MIDPOINTS, N_TEMPERATURES); 
    HeatFluxProfile = zeros(N_MIDPOINTS, 1);
    ConvectionCoef = zeros(N_MIDPOINTS, 1);
    TotalTemperature = zeros(N_MIDPOINTS, 1);

    for i = 2:Mesh.nCurveNodes

        local_Pr = prandtl(i-1);
        local_gamma = gamma(i-1);
        local_cp = cp(i-1);
        local_mu = mu(i-1);
        local_cstar = cstar(i-1);
        local_T_cea = T_cea(i-1);
        local_Ax = Ax(i-1);
        local_Mach = M(i-1);
        local_T_cold = T_ColdWall(i-1);

        T_tot = local_T_cea * (1 + (local_gamma-1)/2 * local_Mach^2);
        TotalTemperature(i-1) = T_tot;

        T_aw_wall = wallThermo.adiabaticWallTemp(local_gamma, local_Pr, local_Mach, T_tot);
        TemperatureProfile(i-1, 1) = T_aw_wall;

        h_conv = wallThermo.convHeatCoefBartz(local_mu, local_cp, Rc_t, local_cstar, local_Pr, grain.design.p0, At,...
                                              local_Ax,T_aw_wall, T_tot, local_gamma, local_Mach);
        ConvectionCoef(i-1) = h_conv;

        R_cond_total = 0;
        for k = 1:Mesh.nLayers
            R_cond_total = R_cond_total + (Mesh.layerThickness(k) / k_layers(k));
        end
        R_conv = 1 / h_conv;
        R_total = R_conv + R_cond_total;

        panelHeatFlux = (T_aw_wall - local_T_cold) / R_total;
        HeatFluxProfile(i-1) = panelHeatFlux;

        % solve for prescribed thermal flux from cold to hot
        T_layers = zeros(1, Mesh.nLayers + 1);
        T_layers(1) = local_T_cold; 
        
        for k = 1:Mesh.nLayers % out to in 
            t = Mesh.layerThickness(k);
            k_coef = k_layers(k);
            deltaT = panelHeatFlux * (t / k_coef);
            
            T_layers(k+1) = T_layers(k) + deltaT; 
        end

        % --- Pack Arrays strictly from Inner to Outer ---
        TemperatureProfile(i-1, 1) = T_aw_wall;
        
        % fliplr reverses T_layers so it becomes: [T_hot, T_interfaces..., T_cold]
        T_layers_InnerToOuter = fliplr(T_layers);
        TemperatureProfile(i-1, 2:end) = T_layers_InnerToOuter; 

        Panel.Temp = T_layers_InnerToOuter;
        Panel.HeatFlux = panelHeatFlux;
        Panel.h_conv = h_conv;
        Panel.midPoint = [x_midPanel(i-1), y_midPanel(i-1)];
        Panels = [Panels; Panel]; %#ok<AGROW>
    end

    Solution.Mesh = Mesh;
    Solution.Panels = Panels;
    Solution.TemperatureProfile = TemperatureProfile;
    Solution.HeatFluxProfile = HeatFluxProfile;
    Solution.ConvectionCoef = ConvectionCoef;
    Solution.MidPointCoords = MidPointCoords;
    Solution.Mach = M;
    Solution.TotalTemperature = TotalTemperature;
    Solution.cp = cp;
    Solution.gamma = gamma;
    Solution.prandtl = prandtl;
    Solution.mu = mu;
    Solution.cstar = cstar;
    Solution.T_cea = T_cea;
    Solution.AxAt = AxAt;

end