function newMesh = setLayers(Mesh, innerLayerThickness , outerLayerThickness)
    newMesh = Mesh; 
    curveNodes = Mesh.curveNodes; 
    curveNormals = Mesh.curveNormals; 
    nNodes = Mesh.nCurveNodes; 
    innerLayerThickness = innerLayerThickness(:)'; 
    outerLayerThickness = outerLayerThickness(:)'; 
    nInnerLayers = size(innerLayerThickness, 2); 
    nOuterLayers = size(outerLayerThickness, 2); 
     
    % Inner Layer addition (along normals pointing inward)
    subStrate = curveNodes; 
    innerLayerNodes = zeros(2, nNodes, nInnerLayers); % same dim as curveNodes 
    for i = 1:nInnerLayers
        displ = layerThickness(i) * curveNormals; 
        subStrate = subStrate + displ; 
        innerLayerNodes(:, :, i) = subStrate; % growing index i moves from outer layer to inner layer 
    end 


    subStrate = curveNodes; 
    outerLayerNodes = zeros(2, nNodes, nOuterLayers); % same dim as curveNodes 
    for i = 1:nOuterLayers
        displ = layerThickness(i) * curveNormals* (-1); 
        subStrate = subStrate + displ; 
        outerLayerNodes(:, :, i) = subStrate; % growing index i moves from outer layer to inner layer 
    end 

    
    % Outer Layer addition (along normals pointing outward)
    
    nLayers = size(layerNodes, 3); 
    newMesh.layerNodes = layerNodes; 
    newMesh.layerThickness = layerThickness; 
    newMesh.nLayers = nLayers; 
 
    

end 