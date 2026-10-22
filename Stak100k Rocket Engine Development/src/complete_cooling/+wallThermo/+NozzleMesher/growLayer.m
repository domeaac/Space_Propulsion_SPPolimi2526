function newMesh = growLayer(Mesh,nNewLayers, newLayerThickness)


    newMesh = Mesh; 

    normals = Mesh.curveNormals;
    currentLayersNodes = Mesh.layerNodes; 
    curveNodes = Mesh.curveNodes; 
    nCurrentLayers = size(Mesh.layerNodes, 3); 

    addLayers = nNewLayers > 0;

    if addLayers
        % grow layers of specified thickness , starting from the last layer

        % specify initial substrate 
        if nCurrentLayers == 0
            substrate = curveNodes; 
        else 
            substrate = currentLayersNodes(:, :, end);
        end 

        % add layers 
        for i=1:nNewLayers
            newLayer = substrate + newLayerThickness(i)*normals; 
            substrate = newLayer; 
            currentLayersNodes(:,:,end+i) = newLayer; 
        end 

    else %remove layer and ignore thickness 

        currentLayersNodes = currentLayersNodes(:,:,1:end-nNewLayers); 
    end 

    nCurrentLayers=size(nCurrentLayers,3); 
    newMesh.layerNodes = currentLayersNodes; 
    newMesh.curveNodes = curveNodes; 
    newMesh.curveNormals = normals; 
    newMesh.nLayers = nCurrentLayers; 


end 