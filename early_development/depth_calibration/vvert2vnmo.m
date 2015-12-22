function [vnmo] = vvert2vnmo(velocity,delta)

vnmo = velocity .* (1 + 2*delta).^0.5;

end
