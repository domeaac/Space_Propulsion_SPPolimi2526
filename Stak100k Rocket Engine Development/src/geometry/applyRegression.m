function [grain, flag] = applyRegression(grain, timestep, rb)

k = timestep*rb/1e3;
grain.x = grain.x + [k; k; -k; -k];
grain.y = grain.y + [0; -k; -k; 0];

if grain.y(2) <= 0 || grain.y(3) <= 0 || grain.x(3) <= grain.x(2)
    flag = true;
else
    flag = false;

end

