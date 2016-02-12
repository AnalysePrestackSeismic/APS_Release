function [AI,VpVsR,Mu,K] = RockProps(Vp,Vs,Rho)
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

