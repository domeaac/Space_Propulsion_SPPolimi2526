function [grain] = initGrainMesh(grain)

grain.x = [0; 0; grain.L0; grain.L0];
grain.y = [0; grain.Re0-grain.Ri0; grain.Re0-grain.Ri0; 0];


end

