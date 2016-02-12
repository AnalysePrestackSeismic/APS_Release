function [Aiso,Biso,Ciso,Aaniso,Baniso,Caniso,RPPshuey_iso,RPPshuey_aniso] = Shuey(MaxAngle,Vp1,Vs1,Rho1,Vp2,Vs2,Rho2,delta1,delta2,epsilon1,epsilon2)
%SHUEY Summary of this function goes here
%   Detailed explanation goes here

Ddelta=delta2-delta1;
Deps=epsilon2-epsilon1;
Mdelta=(delta1+delta2)/2;
Meps=(epsilon1+epsilon2)/2;

%MaxAngle=60;
n=length(Vp1);

% Anisotropic modelling uses the Chan et al correction to the Ruger
% approximation

% Coefficients below are appripriate for Vp/Vs ratio of 1.8
b1=-0.74;
b2=1.04;
c0=1;
c1=8.42;
c2=26.69;
c3=0.86;

for i = 1:n;
    for deg = 1:MaxAngle+1;
        
        Theta1 = (deg-1)*(pi/180);

        Aaniso(i) = 1/2*(((Vp2(i)-Vp1(i))/((Vp2(i)+Vp1(i))/2))+((Rho2(i)-Rho1(i))/((Rho2(i)+Rho1(i))/2)));
        Baniso(i) = ((1/2*((Vp2(i)-Vp1(i))/((Vp2(i)+Vp1(i))/2)))-((2*(((Vs2(i)+Vs1(i))/2)^2)/(((Vp2(i)+Vp1(i))/2)^2))*(((Rho2(i)-Rho1(i))/((Rho2(i)+Rho1(i))/2))+(2*((Vs2(i)-Vs1(i))/((Vs2(i)+Vs1(i))/2))))))+(1/2*Ddelta)+(Aaniso(i)*((b1*(sqrt(Mdelta)))+(b2*(Ddelta))));
        Caniso(i) = (1/2*(((Vp2(i)-Vp1(i))/((Vp2(i)+Vp1(i))/2))+Deps))+((c0*Aaniso(i)^2)+(c1*Aaniso(i)*Meps)+(c2*(Aaniso(i)^2)*Meps)+(c3*(Deps)^2));

        Aiso(i) = 1/2*(((Vp2(i)-Vp1(i))/((Vp2(i)+Vp1(i))/2))+((Rho2(i)-Rho1(i))/((Rho2(i)+Rho1(i))/2)));
        Biso(i) = (1/2*((Vp2(i)-Vp1(i))/((Vp2(i)+Vp1(i))/2)))-((2*(((Vs2(i)+Vs1(i))/2)^2)/(((Vp2(i)+Vp1(i))/2)^2))*(((Rho2(i)-Rho1(i))/((Rho2(i)+Rho1(i))/2))+(2*((Vs2(i)-Vs1(i))/((Vs2(i)+Vs1(i))/2)))));
        Ciso(i) = 1/2*((Vp2(i)-Vp1(i))/((Vp2(i)+Vp1(i))/2));


        Diso((deg-1)+1,1)=Aiso(i)+(Biso(i)*sin(Theta1)^2)+(Ciso(i)*sin(Theta1)^2*tan(Theta1)^2);
        Daniso((deg-1)+1,1)=Aaniso(i)+(Baniso(i)*sin(Theta1)^2)+(Caniso(i)*sin(Theta1)^2*tan(Theta1)^2);
    
    end     
    
        RPPshuey_iso(1:MaxAngle+1,i)=Diso(1:MaxAngle+1,1);
        RPPshuey_aniso(1:MaxAngle+1,i)=Daniso(1:MaxAngle+1,1);
end

