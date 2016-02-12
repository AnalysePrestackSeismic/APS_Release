function [AI,VpVsR,Mu,K] = RockProps(Vp,Vs,Rho)
%ROCKPROPS Summary of this function goes here
%   Detailed explanation goes here

% Bulk rock properties 
% Acoustic Impedance
AI=Vp.*Rho;

% Vp/Vs
VpVsR=Vp./Vs;

% Shear Moduli
Mu=Vs.^2.*Rho;

% Bulk Moduli
K=(Vp.^2.*Rho)-((4/3)*Mu);

end

