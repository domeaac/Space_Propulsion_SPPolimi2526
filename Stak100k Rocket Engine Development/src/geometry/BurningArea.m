function [Area] = BurningArea(grain)

Area = 2*pi*(grain.Re0^2 - (grain.Re0-grain.y(2))^2) + (grain.x(3)-grain.x(2))*2*pi*(grain.Re0-grain.y(2));

end

