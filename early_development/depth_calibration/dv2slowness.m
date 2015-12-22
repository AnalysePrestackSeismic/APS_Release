function [ t_slowness ] = dv2slowness( depth_vint )
%DV2SLOWNESS Convert depth-vint to time-slowness
%   depth_vint: array of time-depth pairs
%   t_slowness: array of time-slowness pairs

sl = depth_vint(:,2).^-1;
dd = depth_vint(2:end,1)-depth_vint(1:end-1,1);
dt = dd.*sl(2:end);
t = [0;cumsum(dt)];

t_slowness = [t,sl];

end

