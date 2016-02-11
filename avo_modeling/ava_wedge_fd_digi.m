function [output] =ava_wedge_fd_digi(trend_file_overburden,trend_file_reservoir,job_meta_path,wavelet_file,depth,method,plot_results)
%% ------------------ Disclaimer  ------------------
% 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) makes no representation or warranty, express or implied, in 
% respect to the quality, accuracy or usefulness of this repository. The code
% is this repository is supplied with the explicit understanding and 
% agreement of recipient that any action taken or expenditure made by 
% recipient based on its examination, evaluation, interpretation or use is 
% at its own risk and responsibility.
% 
% No representation or warranty, express or implied, is or will be made in 
% relation to the accuracy or completeness of the information in this 
% repository and no responsibility or liability is or will be accepted by 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) in relation to it.
%% ------------------ License  ------------------ 
% GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
%% github
% https://github.com/AnalysePrestackSeismic/
%% ------------------ FUNCTION DEFINITION ---------------------------------
[vp1,vs1,rho1]=access_depth_trends(trend_file_overburden,depth);
[vp2,vs2,rho2]=access_depth_trends(trend_file_reservoir,depth);
[output] =ava_wedge_digi(vp1,vs1,rho1,vp2,vs2,rho2,job_meta_path,wavelet_file,depth,method,plot_results);
end
