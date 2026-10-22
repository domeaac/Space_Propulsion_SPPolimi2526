function newMesh = newSetLayers(Mesh, innerLayerThickness, outerLayerThickness)
% SETLAYERS Generates a multi-layered mesh and concatenates thickness arrays.
%
% Inputs:
%   Mesh                - Struct containing curveNodes (2xN), curveNormals (2xN).
%   innerLayerThickness - Row vector of thicknesses for inward layers (+n direction).
%   outerLayerThickness - Row vector of thicknesses for outward layers (-n direction).
%
% Outputs:
%   newMesh             - Updated struct with layerNodes and a unified layerThickness 
%                         array ordered from outermost to innermost.

    newMesh = Mesh; 
    curveNodes = Mesh.curveNodes; 
    curveNormals = Mesh.curveNormals; 
    nNodes = Mesh.nCurveNodes; 
    
    innerLayerThickness = innerLayerThickness(:)'; 
    outerLayerThickness = outerLayerThickness(:)'; 
    
    nInnerLayers = length(innerLayerThickness); 
    nOuterLayers = length(outerLayerThickness); 
     
    subStrateInner = curveNodes; 
    innerLayerNodes = zeros(2, nNodes, nInnerLayers); 
    for i = 1:nInnerLayers
        displ = innerLayerThickness(i) * curveNormals; 
        subStrateInner = subStrateInner + displ; 
        innerLayerNodes(:, :, i) = subStrateInner; 
    end 

    subStrateOuter = curveNodes; 
    outerLayerNodes = zeros(2, nNodes, nOuterLayers); 
    for i = 1:nOuterLayers
        displ = outerLayerThickness(i) * curveNormals * (-1); 
        subStrateOuter = subStrateOuter + displ; 
        outerLayerNodes(:, :, i) = subStrateOuter; 
    end 

    TotNewLayersNodes = cat(3, flip(outerLayerNodes, 3), curveNodes, innerLayerNodes); 
    
    newMesh.layerNodes = TotNewLayersNodes; 
    newMesh.nLayers = size(TotNewLayersNodes, 3) - 1;
    newMesh.layerThickness = [fliplr(outerLayerThickness), innerLayerThickness];
end