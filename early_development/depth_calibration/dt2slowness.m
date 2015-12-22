function [ t_slowness ] = dt2slowness( depth_time )
%TD2SLOWNESS Convert depth-time to time-slowness
%   t_depth: array of time-depth pairs
%   t_slowness: array of time-slowness pairs

if depth_time(1,1) ~= 0
    depth_time = [0,0;depth_time];
end

dt = depth_time(2:end,2)-depth_time(1:end-1,2);
dd = depth_time(2:end,1)-depth_time(1:end-1,1);

t_slowness = [depth_time(2:end,2),dt./dd];

end

