function [ t_depth ] = dv2tdepth( depth_vint )
%DV2SLOWNESS Convert depth-vint to time-depth
%   depth_vint: array of time-depth pairs
%   t-depth: array of time-depth pairs

sl = depth_vint(:,2).^-1;
dd = depth_vint(2:end,1)-depth_vint(1:end-1,1);
dt = dd.*sl(2:end);
t = [0;cumsum(dt)];

t_depth = [t,depth_vint(:,1)];

end

