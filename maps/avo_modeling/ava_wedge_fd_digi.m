function [output] =ava_wedge_fd_digi(trend_file_overburden,trend_file_reservoir,job_meta_path,wavelet_file,depth,method,plot_results)
[vp1,vs1,rho1]=access_depth_trends(trend_file_overburden,depth);
[vp2,vs2,rho2]=access_depth_trends(trend_file_reservoir,depth);
[output] =ava_wedge_digi(vp1,vs1,rho1,vp2,vs2,rho2,job_meta_path,wavelet_file,depth,method,plot_results);
end
