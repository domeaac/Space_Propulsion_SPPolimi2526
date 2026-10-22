function [Volume] = ccVolume(grain)
Volume = (grain.x(3)-grain.x(2))*pi*(grain.Re0-grain.y(2))^2 + pi*(grain.Re0^2) *(2*grain.x(1));
end

