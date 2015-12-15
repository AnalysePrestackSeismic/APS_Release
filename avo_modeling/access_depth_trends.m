function[vp,vs,rho]=access_depth_trends(trend_file,depth)
trend = read_las_file(trend_file);
dz=trend.step;
index=floor(depth/dz)+1;
vp=trend.curves(index,2);
vs=trend.curves(index,3);
rho=trend.curves(index,3);
end
