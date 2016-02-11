function [intercept_reconstruct,gradient_reconstruct,intercept_reconstruct_s,gradient_reconstruct_s] = ava_hodogram(intercept,gradient,freq_c,stretch_factor,phase_rad,plotfig)
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
%% function for plotting ava hodograms
%   INPUTS:
%       intercept: intercept value;
%       gradient: gradient value
%       freq_c = zero offset stack central frequency
%       stretch factor: central frequecy(near stack)/central frequency (far stack)
%       phase_rad : phase rotation of wavelet in radians
%       plotfig : plot figures

%   OUTPUTS:
%       reconstructed intercept and gradient with and without stretching

%   PLOTS:
%       fig1: wavelet at zero offset
%       fig2: ava bhehaviour plots without stretching
%       fig3: ava behaviour with stretching
%       fig4: intercept and gradient traces with and without tracing

% EXAMPLE:
%    ava_hodogram(100,-500,30,2,0,1);

%-------------------------------------------------------------------------------
%%
close all;
n_fold=41;
phase = phase_rad;
angle=0:(n_fold-1);
angle_t=(sind(angle).*sind(angle));
[wavelet,t]=ricker(freq_c,n_s_h,0.001,0,phase);% create wavelet
if plotfig==1
    figure(1);
    plot(t,wavelet);xlabel('time');title('wavelet');
end
r=intercept +angle_t .*gradient;    % create refletivity across angles

% ----------------if wavelet is stationary across offset-------------------------
r_synthetic=conv2(wavelet',r);  %create synthetic gather

intercept_orig=zeros(2*n_s_h-1,1);intercept_orig(n_s_h)=intercept;
gradient_orig=zeros(2*n_s_h-1,1);gradient_orig(n_s_h)=gradient;


intercept_reconstruct=zeros(1,size(r_synthetic,1));
gradient_reconstruct=zeros(1,size(r_synthetic,1));
for n_t=1:size(r_synthetic,1)
    reflections=r_synthetic(n_t,:);
    p=polyfit(angle_t,reflections,1);
    intercept_reconstruct(n_t)=p(2);
    gradient_reconstruct(n_t)=p(1);
end
if plotfig==1
    figure(2);
    subplot(2,2,1);scatter(angle,r);xlabel('angle (degrees)');ylabel('reflectivity');title('AVA of reflectivity');
    subplot(2,2,3);imagesc(r_synthetic);xlabel('angle (degrees)');ylabel('time_sample');title('Synthetic AVA gather')
    subplot(2,2,2);plot(r_synthetic');
    subplot(2,2,4);scatter(intercept_reconstruct,gradient_reconstruct);xlabel('intercept'),ylabel('gradient');title('hodogram with no wavelet stretch')
end



% -----------if wavelet stretches across offset------------------
freq_c_s=freq_c./(1:((stretch_factor-1)/(n_fold-1)):stretch_factor);
r_synthetic_s=zeros((2*n_s_h-1),n_fold);
for i=1:n_fold
    [wavelet_s,t]=ricker(freq_c_s(i),n_s_h,0.001,0,phase);
    r_synthetic_s(:,i)=conv(wavelet_s,r(i));
end

intercept_reconstruct_s=zeros(1,size(r_synthetic_s,1));
gradient_reconstruct_s=zeros(1,size(r_synthetic_s,1));
for n_t=1:size(r_synthetic_s,1)
    reflections=r_synthetic_s(n_t,:);
    p=polyfit(angle_t,reflections,1);
    intercept_reconstruct_s(n_t)=p(2);
    gradient_reconstruct_s(n_t)=p(1);
end
if plotfig==1
    figure(3);
    subplot(2,2,1);scatter(angle,r);xlabel('angle (degrees)');ylabel('reflectivity');title('AVA of reflectivity');
    subplot(2,2,3);imagesc(r_synthetic_s);xlabel('angle (degrees)');ylabel('time_sample');title('Synthetic AVA gather')
    subplot(2,2,2);plot(r_synthetic_s');
    subplot(2,2,4);scatter(intercept_reconstruct_s,gradient_reconstruct_s);xlabel('intercept'),ylabel('gradient');title('hodogram with wavelet stretch')
end

%--------- plot diferent versions of intercept and gradient traces ----------

if plotfig==1
    figure(4);
    subplot(2,3,1); plot(intercept_orig,t,'-r');title ('Original Intercept');xlim([-2*abs(intercept),2*abs(intercept)]);grid('on');
    subplot(2,3,2);plot(intercept_reconstruct,t,'-r');title ('Intercept without stretch');xlim([-2*abs(intercept),2*abs(intercept)]);grid('on');
    subplot(2,3,5);plot(gradient_reconstruct,t);title ('Gradient without stretch');xlim([-2*abs(gradient),2*abs(gradient)]);grid('on');
    subplot(2,3,4); plot(gradient_orig,t);title('Original Gradient');xlim([-2*abs(gradient),2*abs(gradient)]);grid('on');
    subplot(2,3,3);plot(intercept_reconstruct_s,t,'-r');title ('Intercept with stretch');xlim([-2*abs(intercept),2*abs(intercept)]);grid('on');
    subplot(2,3,6);plot(gradient_reconstruct_s,t);title('Gradient with stretch');xlim([-2*abs(gradient),2*abs(gradient)]);grid('on');
end

end
function [s,t] = ricker(f,n,dt,t0,phase)
%RICKER creates an causal ricker wavelet signal
T = dt*(n-1);
t = (-1*T):dt:T;
tau = t-t0;
%s = (1-tau.*tau*f^2*pi^2).*exp(-tau.^2*pi^2*f^2+phase);
s = (1-tau.*tau*f^2*pi^2).*exp(-(tau*pi*f-phase).^2);
end