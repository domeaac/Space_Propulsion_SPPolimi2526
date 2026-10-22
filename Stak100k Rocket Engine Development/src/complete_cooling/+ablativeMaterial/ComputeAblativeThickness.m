function [ablativeThickness, initialThickness, virginThickness] = ComputeAblativeThickness(burnTime, rateOfBurn, alpha_perp, SafetyFactor)
% COMPUTEABLATIVETHICKNESS Computes the required thermal protection system thickness.
%
% Evaluates the total required thickness of an ablative composite wall by 
% summing the material thickness lost to surface erosion and the virgin 
% material depth required for insulation, driven by thermal penetration.
%
% Inputs:
%   burnTime   - Total duration of the thermal exposure [s]
%   rateOfBurn - Material erosion rate across the ply [mm/s]
%   alpha_perp - Thermal diffusivity perpendicular to the ply [cm^2/s]
%
% Outputs:
%   ablativeThickness - Total required TPS coating thickness [mm]
%   initialThickness  - Expected thickness lost to active ablation [mm]
%   virginThickness   - Required insulation thickness (thermal penetration depth) [mm]

    initialThickness = rateOfBurn * burnTime; 
    
    virginThickness = sqrt(alpha_perp * 1e-4 * burnTime); 
    virginThickness = virginThickness * 1000; 
    
    ablativeThickness = initialThickness + virginThickness; 

    ablativeThickness = ablativeThickness * SafetyFactor; 





    fprintf('\n--- Ablative TPS Sizing Log ---\n');
    fprintf('Burn Time:                 %.2f s\n', burnTime);
    fprintf('Erosion Rate (v_er):       %.3f mm/s\n', rateOfBurn);
    fprintf('Thermal Diffusivity (a):   %.4f cm^2/s\n', alpha_perp);
    fprintf('----------------------------------------\n');
    fprintf('Eroded Thickness:          %.3f mm\n', initialThickness);
    fprintf('Virgin Insulation:         %.3f mm\n', virginThickness);
    fprintf('Safety Factor:             %.1f \n', SafetyFactor);
    fprintf('Total Required Thickness:  %.3f mm\n', ablativeThickness);
    fprintf('----------------------------------------\n\n');

end