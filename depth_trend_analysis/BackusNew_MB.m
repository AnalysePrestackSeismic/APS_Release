%Backus average
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


%IMPORTANT:   a should not be zero

function [Vpn,Vsn,RHOBn]=BackusNew_MB(Vp,Vs,RHOB,n)


RHOBn=ArithmNew_MB(RHOB,n,1);
M=Vp.*Vp.*RHOB;
G=Vs.*Vs.*RHOB;

Mn=HarmonicNew_MB(M,n);
Gn=HarmonicNew_MB(G,n);

PRn=.5*(Mn./Gn-2)./(Mn./Gn-1);
Ipn=sqrt(Mn.*RHOBn);
Isn=sqrt(Gn.*RHOBn);

Vpn=sqrt(Mn./RHOBn);
Vsn=sqrt(Gn./RHOBn);

