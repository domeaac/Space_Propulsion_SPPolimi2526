function croppedMesh = cropMesh(Mesh, AxAtSpan, cutType)
% CROPMESH Extracts a specific axial segment of a nozzle mesh based on Area Ratio (A/At).
%
% Inputs:
%   Mesh     - (1,1) struct containing the full nozzle mesh arrays.
%   AxAtSpan - (1,2) double array specifying the start and end A/At values.
%   cutType  - (1,1) string specifying the nozzle region to crop:
%              "sub"     : Subsonic (convergent) region
%              "sup"     : Supersonic (divergent) region
%              "sub-sup" : Transonic (convergent-divergent) spanning the throat

    arguments
        Mesh (1,1) struct
        AxAtSpan (1,2) double
        cutType (1,1) string {mustBeMember(cutType, ["sub", "sup", "sub-sup"])}
    end

    ARatio1 = AxAtSpan(1);
    ARatio2 = AxAtSpan(2);

    croppedMesh = Mesh;
    curveNodes = Mesh.curveNodes;
    layerNodes = Mesh.layerNodes;
    curveNormals = Mesh.curveNormals;
    At = Mesh.At;

    % Compute nodal area ratios
    y_nodes = curveNodes(2,:);
    A_nodes = pi * y_nodes.^2;
    AxAtNodes = A_nodes ./ At;

    min_AcutRatio = min(AxAtSpan); 
    max_AcutRatio = max(AxAtSpan); 
    maxAreaRatio = max(AxAtNodes);
    minAreaRatio = min(AxAtNodes); % Expected to be ~1.0 at the throat

    % Corrected boundary logic
    if min_AcutRatio < minAreaRatio || max_AcutRatio > maxAreaRatio
        error("cropMesh:OutOfBounds", ...
              "Requested cut span [%.4g, %.4g] is out of bounds. Available ratio range is [%.4g, %.4g].", ...
              min_AcutRatio, max_AcutRatio, minAreaRatio, maxAreaRatio);
    end 

    % Locate the aerodynamic throat to split the monotonic domains
    [~, throatID] = min(AxAtNodes);

    switch cutType
        case 'sub'
            % Convergent section: Ratio must be decreasing
            if ARatio1 < ARatio2
                error("cropMesh:IllegalSpan", "Illegal cut: ARatio1 < ARatio2 for a convergent section.");
            end
            
            % Restrict search strictly to the convergent side (1 to throat)
            firstCutID = find(AxAtNodes(1:throatID) <= ARatio1, 1, "first");
            lastCutID  = find(AxAtNodes(1:throatID) <= ARatio2, 1, "first");

        case 'sup'
            % Divergent section: Ratio must be increasing
            if ARatio1 > ARatio2
                error("cropMesh:IllegalSpan", "Illegal cut: ARatio1 > ARatio2 for a divergent section.");
            end
            
            % Restrict search strictly to the divergent side (throat to end)
            firstCutID = throatID - 1 + find(AxAtNodes(throatID:end) >= ARatio1, 1, "first");
            lastCutID  = throatID - 1 + find(AxAtNodes(throatID:end) >= ARatio2, 1, "first");

        case 'sub-sup'
            % Transonic section: ARatio1 in convergent, ARatio2 in divergent
            firstCutID = find(AxAtNodes(1:throatID) <= ARatio1, 1, "first");
            lastCutID  = throatID - 1 + find(AxAtNodes(throatID:end) >= ARatio2, 1, "first");
            
    end

    if isempty(firstCutID) || isempty(lastCutID) || (firstCutID > lastCutID)
        error("cropMesh:NodesNotFound", "Requested area ratio span could not be located in the mesh boundaries.");
    end

    % Slice geometry matrices (Note: layerNodes is sliced along 2nd dimension)
    curveNodes = curveNodes(:, firstCutID:lastCutID);
    layerNodes = layerNodes(:, firstCutID:lastCutID, :); 
    curveNormals = curveNormals(:, firstCutID:lastCutID);
    nCurveNodes = size(curveNodes, 2);
    
    % Write to output struct
    croppedMesh.curveNodes = curveNodes;
    croppedMesh.layerNodes = layerNodes;
    croppedMesh.curveNormals = curveNormals;
    croppedMesh.nCurveNodes = nCurveNodes;
    croppedMesh.region = cutType;

end