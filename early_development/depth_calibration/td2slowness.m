function [ t_slowness ] = td2slowness( t_depth )
%TD2SLOWNESS Convert time-depth to time-slowness
%   t_depth: array of time-depth pairs
%   t_slowness: array of time-slowness pairs

if t_depth(1,1) ~= 0
    t_depth = [0,0;t_depth];
end

dt = t_depth(2:end,1)-t_depth(1:end-1,1);
dd = t_depth(2:end,2)-t_depth(1:end-1,2);

t_slowness = t_depth(2:end,dt./dd);

end

