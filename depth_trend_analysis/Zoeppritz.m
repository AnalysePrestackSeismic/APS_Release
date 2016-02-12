function [RPP,RPS,TPP,TPS] = Zoeppritz(MaxAngle,Vp1,Vs1,Rho1,Vp2,Vs2,Rho2);
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
%ZOEPPRITZ Summary of this function goes here
%   Detailed explanation goes here

%MaxAngle=60;
l2=MaxAngle*4;
n=length(Vp1);

for i = 1:n
    for deg = 1:MaxAngle+1
    Theta1 = (deg-1)*(pi/180);
    Theta2 = asin(Vp2(i)*(sin(Theta1)/Vp1(i)));
    Sigma1 = asin(Vs1(i)*(sin(Theta1)/Vp1(i)));
    Sigma2 = asin(Vs2(i)*(sin(Theta1)/Vp1(i)));
    
  
   c1=2*Rho1(i)*Vs1(i);
   c2=Rho1(i)*Vs1(i);
   c3=2*Rho2(i)*Vs2(i);
   c4=Rho2(i)*Vs2(i);
   d1=Rho1(i)*Vp1(i);
   d2=Rho1(i)*Vs1(i);
   d3=Rho2(i)*Vp2(i);
   d4=Rho2(i)*Vs2(i);
  
   
   G= [-sin(Theta1) -cos(Sigma1) sin(Theta2) cos(Sigma2);
       cos(Theta1) -sin(Sigma1) cos(Theta2) -sin(Sigma2);
       c1*sin(Sigma1)*cos(Theta1) c2*(1-2*sin(Sigma1)^2) c3*sin(Sigma2)*cos(Theta2) c4*(1-2*sin(Sigma2)^2);
       -d1*(1-2*sin(Sigma1)^2) d2*sin(2*Sigma1) d3*(1-2*(sin(Sigma2)^2)) -d4*sin(2*Sigma2)];
   
   D= [sin(Theta1) cos(Sigma1) -sin(Theta2) -cos(Sigma2);
       cos(Theta1) -sin(Sigma1) cos(Theta2) -sin(Sigma2);
       c1*sin(Sigma1)*cos(Theta1) c2*(1-2*sin(Sigma1)^2) c3*sin(Sigma2)*cos(Theta2) c4*(1-2*sin(Sigma2)^2);
       d1*(1-2*sin(Sigma1)^2) -d2*sin(2*Sigma1) -d3*(1-2*(sin(Sigma2)^2)) d4*sin(2*Sigma2)];
   
   
    M(4*(deg-1)+1:4*(deg-1)+4,1:4)=G\D;
    %RPPa(1:deg,1)=((delta/2)*sin(Theta1)^2)+((epsilon/2)*sin(Theta1)^2*tan(Theta1)^2);
    end
    
    RPP(1:MaxAngle+1,i)=M(1:4:l2+1,1);
    RPS(1:MaxAngle+1,i)=M(2:4:l2+2,1);
    TPP(1:MaxAngle+1,i)=M(3:4:l2+3,1);
    TPS(1:MaxAngle+1,i)=M(4:4:l2+4,1);
    %RPPanis=RPP+RPPa;
end
   
end

