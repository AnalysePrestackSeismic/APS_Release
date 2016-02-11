function [x,theta,twt,z_dec] = snellr_tracer(file_path)
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

% INPUT: file path : path to exported velocity file 

% OUTPUT: 
%   x: offset reached
%   theta: angle of entry
%   twt : two way time
%   z_dec : z axis

%--------------------------------------------------------------
% read the log exported
log_read=dlmread(file_path);
figure(1);
plot(log_read(:,2),log_read(:,1));set(gca,'YDir','rev');
hold on;

z_max=max(log_read(:,1));
z_log=log_read(:,1);
vel_log=log_read(:,2);

% decimate the velocity profile and do averaging within decimation
s_dec=200;
n_dec=1:s_dec:length(z_log);
z_dec=z_log(n_dec);
n_max=length(z_dec);
vel_dec=zeros(n_max,1);
for i1=1:(n_max-1)
vel_dec(i1)=rms(vel_log(n_dec(i1):(n_dec(i1+1)-1)));
end
vel_dec(n_max)=vel_dec(n_max-1);
z_s=z_dec(2)-z_dec(1);

scatter(vel_dec,z_dec);set(gca,'YDir','rev');
hold off;

% ray tracing
theta0=10;
x=zeros(n_max,1);
theta=zeros(n_max,1);
time_ow=zeros(n_max,1);
theta(1)=theta0;
x(1)=0;
time_ow(1)=0;

for i2=2:n_max
    theta(i2)=real(asind((vel_dec(i2)/vel_dec(i2-1))*sind(theta(i2-1))));
    dx=z_s*tand(theta(i2));
    x(i2)=x(i2-1)+dx;
    time_ow(i2)=time_ow(i2-1)+(((z_s^2+dx^2)^0.5)/vel_dec(i2));
end
twt=2*time_ow;
offset=2*x;
figure(2);
plot(x,z_dec);set(gca,'YDir','rev');
figure(3);
plot(theta,z_dec);set(gca,'YDir','rev');
figure(4);
plot(twt,z_dec);set(gca,'YDir','rev');

end