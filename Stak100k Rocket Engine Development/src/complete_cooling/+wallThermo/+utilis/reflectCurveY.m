function [nodesFlipped, normalsFlipped] = reflectCurveY(nodes, normals)
% REFLECTCURVEY Reflects a 2D curve and its associated normal vectors across the Y-axis.
%
% Inputs:
%   nodes   : 2xN array of curve coordinates [x; y]
%   normals : 2xN array of normal vectors [nx; ny]
%
% Outputs:
%   nodesFlipped   : 2xN array of reflected and reordered coordinates
%   normalsFlipped : 2xN array of reflected and reordered normals

    % 1. Geometric Reflection (X -> -X)
    nodesFlipped = nodes;
    nodesFlipped(1,:) = -nodesFlipped(1,:);
    
    
    normalsFlipped = normals; %noFlip
    
    
    % 3. Restore Parameterization Direction
    nodesFlipped = fliplr(nodesFlipped);
    normalsFlipped = fliplr(normalsFlipped);
end