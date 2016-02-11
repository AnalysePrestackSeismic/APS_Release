function [angle,Refl] =avo_3term_aki_richard(vp1,vs1,rho1,vp2,vs2,rho2,max_angle)
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
%% FUNCTION TO COMPUTE Half space elastic reflectivity .

plot_on=0;close all;
drho=rho2-rho1;
dvp=vp2-vp1;
dvs=vs2-vs1;

avg_rho=(rho2+rho1)/2;
avg_vp=(vp2+vp1)/2;
avg_vs=(vs2+vs1)/2;


a=0.5*(dvp/avg_vp+drho/avg_rho);
b=0.5*dvp/avg_vp - 4*((avg_vs/avg_vp).^2).*(dvs/avg_vs)-2*((avg_vs/avg_vp).^2).*(drho/avg_rho);
c=0.5*dvp/avg_vp;

angle=0:1:max_angle;
Refl = a+b*((sind(angle)).^2)+c*((tand(angle)).^2).*((sind(angle)).^2);
if plot_on==1
    scatter(angle,Refl);
    xlabel('Angle'); ylabel('Refletivity');
end
end
