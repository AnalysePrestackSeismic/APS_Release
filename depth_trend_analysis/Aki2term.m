function [R0,G] = Aki2term(Vp1,Vs1,Rho1,Vp2,Vs2,Rho2)
%AKI-RICHARDS-2TERM Summary of this function goes here
%   Detailed explanation goes here

n=length(Vp1);

for i = 1:n 

R0(i,1) = 1/2*(((Vp2(i)-Vp1(i))/((Vp2(i)+Vp1(i))/2))+((Rho2(i)-Rho1(i))/((Rho2(i)+Rho1(i))/2)));
G(i,1) = (1/2*((Vp2(i)-Vp1(i))/((Vp2(i)+Vp1(i))/2)))-((2*(((Vs2(i)+Vs1(i))/2)^2)/(((Vp2(i)+Vp1(i))/2)^2))*(((Rho2(i)-Rho1(i))/((Rho2(i)+Rho1(i))/2))+(2*((Vs2(i)-Vs1(i))/((Vs2(i)+Vs1(i))/2)))));

end

