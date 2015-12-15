function [ava_refl_top,time_thickness_true,time_thickness_ap] =ava_wedge_fd(trend_file_overburden,trend_file_reservoir,depth,max_angle,freq_c,method,plot_results)
[vp1,vs1,rho1]=access_depth_trends(trend_file_overburden,depth);
[vp2,vs2,rho2]=access_depth_trends(trend_file_reservoir,depth);
[ava_refl_top,time_thickness_true,time_thickness_ap] =ava_wedge(vp1,vs1,rho1,vp2,vs2,rho2,max_angle,freq_c,method,plot_results);
end
