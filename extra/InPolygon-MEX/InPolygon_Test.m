
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
% Test inpolygon 
close all
figure
[X,Y]=meshgrid([0:0.1:5],[0:0.1:2]);

XV = [0 5 5 3 4 1 2 0 0];
YV = [0 0 2 2 1 1 2 2 0];
tic
[in_on,on,strict] = InPolygon(X,Y,XV,YV);
toc
plot(XV,YV,X(strict),Y(strict),'r+',X(~in_on),Y(~in_on),'b+',X(on),Y(on),'go')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rand('state',0)
dtheta = pi/10;
theta  = (-pi:dtheta:(pi-dtheta))';
node   = [cos(theta) sin(theta)];

% Test points
p = 3*(rand(1e5,2)-0.5);

px=p(:,1);
py=p(:,2);
nx=node(:,1);
ny=node(:,2);
tic, [in_on,on,strict] =  InPolygon(px',py',nx,ny); t3 = toc;
 H Hfigure
plot(nx,ny,px(strict),py(strict),'r+',px(~in_on),py(~in_on),'b+',px(on),py(on),'go')

