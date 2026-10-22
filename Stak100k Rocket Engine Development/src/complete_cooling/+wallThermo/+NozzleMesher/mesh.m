function Mesh = mesh(NozzleCurve, curveDomain, nLayers, layerThickness, ds)
% MESH Generates a structured layered mesh from a 1D base parametric curve.
%
% Inputs:
%   NozzleCurve    - Function handle defining the nozzle profile y = g(x)
%   curveDomain    - 1x2 double array [x_min, x_max] for the spatial domain
%   nLayers        - Integer specifying the number of extruded layers
%   layerThickness - 1D double array of length nLayers specifying thickness
%   ds             - Double specifying the arc-length discretization step
%
% Outputs:
%   Mesh           - Struct containing curveNodes, layerNodes, curveNormals nLayers and layerthicness

    arguments
        NozzleCurve (1,1) function_handle; 
        curveDomain (1,2) double; 
        nLayers (1,1) uint16; 
        layerThickness (1,:) double; 
        ds (1,1) double; 
    end 

    if nLayers ~= length(layerThickness)
        error("mesh:IllegalLayer", "Number of layers must match the layerThickness array length.");
    end 

    nSamples = 10000; 
    x_fine = linspace(curveDomain(1), curveDomain(2), nSamples); 
    y_fine = NozzleCurve(x_fine);
    
    dx = diff(x_fine);
    dy = diff(y_fine); 
    segmentLengths = sqrt(dx.^2 + dy.^2);
    cumLen = [0, cumsum(segmentLengths)]; 
    totLen = cumLen(end); 
    
    targetLen = 0:ds:totLen; 
    if targetLen(end) < totLen
        targetLen(end+1) = totLen; 
    end 

    x_nodes = interp1(cumLen, x_fine, targetLen, 'pchip'); % this works 
    y_nodes = NozzleCurve(x_nodes); 

    curveNodes = [x_nodes(:)'; y_nodes(:)']; % 2xN array 
    nNodes = size(curveNodes, 2); 

    t_hat_seg = diff(curveNodes, 1, 2); %tangent dir 
    t_hat_seg = t_hat_seg ./ vecnorm(t_hat_seg, 2, 1); %normalized
    n_hat_seg = [-t_hat_seg(2, :); t_hat_seg(1, :)]; % flipped to get normal 

    flip_idx = n_hat_seg(2, :) > 0;
    n_hat_seg(:, flip_idx) = -n_hat_seg(:, flip_idx); %flip if normal point UP 

    curveNormals = zeros(2, nNodes);
    curveNormals(:, 1) = n_hat_seg(:, 1);
    curveNormals(:, end) = n_hat_seg(:, end);
    for i = 2:(nNodes-1)
        n_avg = n_hat_seg(:, i-1) + n_hat_seg(:, i); %avg per point normal 
        curveNormals(:, i) = n_avg / norm(n_avg); % normalize 
    end 

    % LAYERS 
    subStrate = curveNodes; 
    layerNodes = zeros(2, nNodes, nLayers); % same dim as curveNodes 
    for i = 1:nLayers
        displ = layerThickness(i) * curveNormals; 
        subStrate = subStrate + displ; 
        layerNodes(:, :, i) = subStrate; 
    end 

    Mesh.layerNodes = layerNodes; 
    Mesh.curveNodes = curveNodes; 
    Mesh.curveNormals = curveNormals; 
    Mesh.nLayers = nLayers; 
    Mesh.nCurveNodes = nNodes;  
    Mesh.layerThickness = layerThickness; 
    Mesh.NozzleCurveHandle = NozzleCurve; 
 

end