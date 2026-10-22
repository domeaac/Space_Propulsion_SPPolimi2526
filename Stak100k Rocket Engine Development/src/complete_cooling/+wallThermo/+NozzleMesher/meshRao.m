function Mesh = meshRao(alpha_div , beta_conv , At, Ae, Rcc, Rao_thickness, N, Lcc , nLayers , layerThickness,options)
    arguments
        alpha_div (1,1) double; 
        beta_conv (1,1) double; 
        At (1,1) double; 
        Ae (1,1) double; 
        Rcc (1,1) double; 
        Rao_thickness (1,1) double; 
        N (1,1) double; 
        Lcc (1,1) double; 
        nLayers (1,1) double; 
        layerThickness (1,:) double; 
        options.PlotNozzle (1,1) logical = false; 
    end 

    if nLayers ~= size(layerThickness(:)',2)
        error("meshRao:IllegalThick", "inconsistent layer specification, nLayers different from nThickness")
    end 


    % Call rao for nodes and normals 
  
    [RaoHandleFunc, r_nodes, ~, N, ~,~,R_curvature_Throat] = wallThermo.RAO.RAO2(alpha_div, beta_conv, Ae, At, Rcc, Rao_thickness, N, Lcc, options.PlotNozzle); 
    
  
    curveNodes = r_nodes; 
    nNodes = size(curveNodes, 2); 
    curveNormals = N; 

    % reflect curve and normals
    %[curveNodes , curveNormals] = wallThermo.utilis.reflectCurveY(curveNodes, curveNormals); 
    % flip norm along Y 
    curveNormals = -curveNormals; 

    % for the curve, find x_throat location 
    [~ , throatID] = min(abs(curveNodes(1,:))); 
    x_throat = curveNodes(1, throatID); 
    y_throat = curveNodes(2, throatID);
    throatPos = [x_throat, y_throat]; 

    % compute Throat area 
    At = pi*y_throat^2; % assuming circ section 

    % Layer addition 
    subStrate = curveNodes; 
    layerNodes = zeros(2, nNodes, nLayers); % same dim as curveNodes 
    for i = 1:nLayers
        displ = layerThickness(i) * curveNormals; 
        subStrate = subStrate + displ; 
        layerNodes(:, :, i) = subStrate; 
    end 

    % x safety 
    nLayers = size(layerNodes, 3); 
   
    Mesh.layerNodes = layerNodes; 
    Mesh.curveNodes = curveNodes; 
    Mesh.curveNormals = curveNormals; 
    Mesh.nLayers = nLayers; 
    Mesh.nCurveNodes = nNodes;  
    Mesh.layerThickness = layerThickness; 
    Mesh.NozzleFuncHandle = RaoHandleFunc; 
    Mesh.CurvatureRadiusThroat = R_curvature_Throat; 
    Mesh.throatPos = throatPos; 

    % bool to check if mesh was cropped 
    Mesh.cropped = false; 
    Mesh.At = At; 
    Mesh.region = "sub-sup"; % means this mesh has both sub and sup parts 
end 