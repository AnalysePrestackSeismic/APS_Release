function [angle,Refl] =avo_zeoppritz(vp1,vs1,rho1,vp2,vs2,rho2,max_angle)
% From The Rock Physics Manual, page 61
%	Written by Jack Dvorkin 06/04

Ip1=vp1.*rho1; Ip2=vp2.*rho2;
PR1=.5*((vp1./vs1).^2-2)./((vp1./vs1).^2-1);
PR2=.5*((vp2./vs2).^2-2)./((vp2./vs2).^2-1);

%	Angle in degrees
tD1=[0:1:max_angle]'; tDZ=tD1;
angle=tD1';
%	Angle in radian
t1=tD1/57.32;

%	Normal reflectivity
R0=(Ip2-Ip1)./(Ip2+Ip1);

%	Reflection and transmission angles
t2=asin(vp2.*sin(t1)./vp1);
p1=asin(vs1.*sin(t1)./vp1);
p2=asin(vs2.*sin(t1)./vp1);

%	Various parameters
a=rho2.*(1-2.*sin(p2).*sin(p2))-rho1.*(1-2.*sin(p1).*sin(p1));
b=rho2.*(1-2.*sin(p2).*sin(p2))+2.*rho1.*sin(p1).*sin(p1);
c=rho1.*(1-2.*sin(p1).*sin(p1))+2.*rho2.*sin(p2).*sin(p2);
d=2.*(rho2.*vs2.*vs2-rho1.*vs1.*vs1);
p=sin(t1)./vp1;
E=b.*cos(t1)./vp1+c.*cos(t2)./vp2;
F=b.*cos(p1)./vs1+c.*cos(p2)./vs2;
G=a-d.*cos(t1).*cos(p2)./(vp1.*vs2);
H=a-d.*cos(t2).*cos(p1)./(vp2.*vs1);
D=E.*F+G.*H.*p.*p;

%	Reflectivity versus angle
RppZ=(F.*(b.*cos(t1)./vp1-c.*cos(t2)./vp2)-H.*p.*p.*(a+d.*cos(t1).*cos(p2)./(vp1.*vs2)))./D;
% RpsZ=(-2.*(cos(t1)./vp1).*(a.*b+c.*d.*cos(t2).*cos(p2)./(vp2.*vs2)).*p.*vp1)./(vs1.*D);

% RppZ=real(RppZ+.00001); RpsZ=real(RpsZ+.00001);
RppZ=real(RppZ);
Refl=RppZ';