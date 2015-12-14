function [angle,Refl] =avo_3term_aki_richard(vp1,vs1,rho1,vp2,vs2,rho2,max_angle)
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
