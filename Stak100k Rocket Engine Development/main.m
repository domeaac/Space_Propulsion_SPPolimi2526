clear; clc; close all; warning('off', 'all')

PathsInit

%% ENGINE DESIGN

engine_graindesignFinal

%% MONTECARLO ANALYSIS 

engine_resultsMontecarlo

%% THERMAL COATING DESIGN + ABLATIVE PROTECTIONS DESIGN + SUPPORT SYSTEM SIZING

engine_thermalCoatingDesign
% be aware: the thermal coating and ablative protection related plots are shown as selectable in one figure ("figures")

% SUPPORT SYSTEM SIZING - DEPENDENT ON THERMAL DESIGN therefore in the same section
hydraulics % this script needs to run after engine_thermalCoatingDesign

