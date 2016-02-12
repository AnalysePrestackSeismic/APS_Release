function [R0,G] = Aki2term(Vp1,Vs1,Rho1,Vp2,Vs2,Rho2)
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
%AKI-RICHARDS-2TERM Summary of this function goes here
%   Detailed explanation goes here

n=length(Vp1);

for i = 1:n 

R0(i,1) = 1/2*(((Vp2(i)-Vp1(i))/((Vp2(i)+Vp1(i))/2))+((Rho2(i)-Rho1(i))/((Rho2(i)+Rho1(i))/2)));
G(i,1) = (1/2*((Vp2(i)-Vp1(i))/((Vp2(i)+Vp1(i))/2)))-((2*(((Vs2(i)+Vs1(i))/2)^2)/(((Vp2(i)+Vp1(i))/2)^2))*(((Rho2(i)-Rho1(i))/((Rho2(i)+Rho1(i))/2))+(2*((Vs2(i)-Vs1(i))/((Vs2(i)+Vs1(i))/2)))));

end

