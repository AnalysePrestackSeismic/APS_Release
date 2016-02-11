function [intercept_reconstruct,gradient_reconstruct,intercept_reconstruct_s,gradient_reconstruct_s] = ava_hodogram_digi(job_meta_path,wavelet_file,st_vol,inc_vol,end_vol,level,intercept,gradient,phase_rad,plotfig)
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
%       job_meta_path = path to job_meta file
%       wavelet_file = name of wavelet file
%       st_vol = # start volume (not angle) , specify as 0 if you want to use the full ange as available per job meta file 
%       inc_vol = # inc volume (not angle) , specify as 0 if you want to use the full ange as available per job meta file
%       end_vol = # end volume (not angle) , specify as 0 if you want to use the full ange as available per job meta file 
%       level = depth in window levels (# of depth window from DIGI)
%       intercept = AVA intercept value
%       gradient = AVA gradient value
%       phase_rad = not currently used = 0
%       plotfig = 1 if you want to plot results, otherwise 0.

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
%% Function to plot ava hodogram based on a given intercept and gradient
freq_c =30;
job_meta = load(job_meta_path);
wavelet_path = strcat(job_meta.wav_directory,wavelet_file);
wavelet_digi=load(wavelet_path);

st_vol=str2double(st_vol);
inc_vol=str2double(inc_vol);
end_vol=str2double(end_vol);
level=str2double(level);
intercept=str2double(intercept);
gradient=str2double(gradient);
phase_rad=str2double(phase_rad);
plotfig=str2double(plotfig);


if st_vol==0 && inc_vol==0 && end_vol==0
    % angles unitialized get from job meta file   
    st_ang=job_meta.tkey_min;
    inc_ang=job_meta.tkey_inc;
    end_ang = job_meta.tkey_max;
   
    st_vol=1;
    inc_vol=1;
    end_vol=1+(end_ang-st_ang)/inc_ang;
else
    %if initialized
   ang= job_meta.tkey_min:job_meta.tkey_inc:job_meta.tkey_max;
   st_ang=ang(st_vol);
%    inc_vol=job_meta.tkey_inc;
   end_ang=ang(end_vol);
   inc_ang =ang(2)-ang(1);
end


%close all;
n_s_h=job_meta.ns_win/2;
phase = phase_rad;
angles=st_ang:inc_ang:end_ang;
angles_plot=angles;
n_fold=length(angles);
angle_t=(sind(angles).*sind(angles));
[wavelet,t]=ricker(freq_c,n_s_h,job_meta.s_rate/1e6,0,phase);% create wavelet
t_fine=min(t):(t(2)-t(1))/10:max(t);% time axis to facilitate interpolation
% if plotfig==1
%     figure(1);
%     plot(t,wavelet);xlabel('time');title('wavelet');
% end
r=intercept +angle_t.*gradient;    % create refletivity across angles

% ----------------if wavelet is stationary across offset-------------------------
r_synthetic=conv2(wavelet',r);  %create synthetic gather

intercept_orig=zeros(2*n_s_h,1);intercept_orig(n_s_h)=intercept;
gradient_orig=zeros(2*n_s_h,1);gradient_orig(n_s_h)=gradient;

[intercept_reconstruct,gradient_reconstruct] =ig_reconstruct_linefit(r_synthetic,angles);

if plotfig==1
    figure(2); 
    subplot(2,2,1);scatter(angles,r);xlabel('angle (degrees)');ylabel('reflectivity');title('AVA of reflectivity');
    subplot(2,2,3);imagesc(r_synthetic);xlabel('angle (degrees)');ylabel('time_sample');title('Synthetic AVA gather')
    subplot(2,2,2);plot(angles,r_synthetic');title('Reconstructed AVA curves')
    subplot(2,2,4);scatter(intercept_reconstruct,gradient_reconstruct);xlabel('intercept'),ylabel('gradient');title('hodogram with no wavelet stretch')
end


% a plot showing the sample values colouring the hodogram
% if plotfig==1
%     figure(11); 
%     %colcj = [0 0 (143/255); (128/255) 1 (128/255); 0 0 (143/255)];
%     %[ cjrgb ] = make_colormap( [0 0.75 2 3 4 5 6 7 8],[0 0 128 255 255 255 128 0 0],[0 0 255  255 0 255 255 0 0],[143 255 128  0 0 0 128 255 143],129,255 );
%     %[ cjrgb ] = make_colormap( [0 1 2  3  4 5 6],[0 0 128 255  128 0 0],[0 0 255   0  255 0 0],[143 255 128   0  128 255 143],129,255 );
%     [ cjrgb ] = make_colormap( [0 1 2  3  4 5 6],[0 0 64 255  64 0 0],[0 0 255   0  255 0 0],[143 255 64   0  64 255 143],127,255 );
%     %colormap(jet)
%     colormap(cjrgb)
%     %subplot(2,2,1);scatter(angles,r);xlabel('angle (degrees)');ylabel('reflectivity');title('AVA of reflectivity');
%     subs = 45;
%     ends = 85;
% 
%     r_syntheticb = r_synthetic(subs:ends,:);
%     intercept_reconstructb = intercept_reconstruct(subs:ends);
%     gradient_reconstructb = gradient_reconstruct(subs:ends);
%     timeb = t(subs:ends);
%     
%     samples = repmat((1:size(r_syntheticb,1))',1,20);
%     
%     
%     %subplot(1,5,2);plot(gradient_reconstructb,timeb);title('Gradient with real wavelet');xlim([-2*abs(gradient),2*abs(gradient)]);grid('on');   
%     xshifts = (150:175);
%     subplot(1,3,1);imagesc(xshifts,timeb,samples);xlabel('');ylabel('time sample');title('sample colour')
%     hold on;
%     plot(intercept_reconstructb,timeb,'-r');title ('Intercept using single ricker wavelet','FontSize',14);xlim([-2*abs(intercept),2*abs(intercept)]);grid('on');
%     hold off;
% 
%     xshifts = (300:350);
%     subplot(1,3,2);imagesc(xshifts,timeb,samples);xlabel('');ylabel('time sample');title('sample colour')
%     hold on;
%     plot(gradient_reconstructb,timeb);title('Gradient using single ricker wavelet','FontSize',14);xlim([-2*abs(gradient),2*abs(gradient)]);grid('on');
%     %set(p,'Color','black','LineWidth',5)
%     hold off;
%     ptsize = 80;
%     cols = linspace(1,10,length(intercept_reconstructb));
%     subplot(1,3,3);scatter(intercept_reconstructb,gradient_reconstructb,ptsize,cols,'fill','s');xlabel('intercept'),ylabel('gradient');title('IG crossplot showing hodogram without wavelet stretch','FontSize',14);grid('on')
% 
% end

% if plotfig==1
%     figure(12);
%     %colormap('bone')
%     subplot(1,5,1);scatter(angles,r);xlabel('angle (degrees)');ylabel('reflectivity');title('synthetic reflectivity');
%     subplot(1,5,2:3);seisplot(r_syntheticb,20,timeb,angles_plot);xlabel('angle (degrees)');ylabel('time sample');title('Synthetic AVA gather using single ricker wavelet');grid('on');
%     xlim([0 ((angles_plot(end)*20)+200)]);
%     %set(gca,'XTick',-pi:pi/2:pi)
%     set(gca,'XTickLabel',st_ang:5:end_ang+20);
%     %subplot(1,4,2);imagesc(r_syntheticb);xlabel('angle (degrees)');ylabel('time sample');title('Synthetic AVA gather using single ricker wavelet')
%     subplot(1,5,4);plot(intercept_reconstructb,timeb,'-r');title ('Intercept from synthetic gather');xlim([-2*abs(intercept),2*abs(intercept)]);grid('on');
%     subplot(1,5,5);plot(gradient_reconstructb,timeb);title('Gradient from synthetic gather');xlim([-2*abs(gradient),2*abs(gradient)]);grid('on');
% end

% a plot showing the sample values colouring the hodogram
% if plotfig==1
%     figure(20); 
%     %subplot(2,2,1);scatter(angles,r);xlabel('angle (degrees)');ylabel('reflectivity');title('AVA of reflectivity');
%     subs = 57;
%     ends = 74;
%     r_syntheticb = r_synthetic(subs:ends,:);
%     intercept_reconstructb = intercept_reconstruct(subs:ends);
%     gradient_reconstructb = gradient_reconstruct(subs:ends);
%     subplot(1,3,1);imagesc(r_syntheticb);xlabel('angle (degrees)');ylabel('time sample');title('Synthetic AVA gather')
%     %samples = repmat((1:size(r_synthetic,1))',1,size(r_synthetic,2));
%     samples = repmat((1:size(r_syntheticb,1))',1,10);
%     subplot(1,3,2);imagesc(samples);xlabel('');ylabel('time sample');title('sample colour')
%     ptsize = 50;
%     cols = linspace(1,10,length(intercept_reconstructb));
%     subplot(1,3,3);scatter(intercept_reconstructb,gradient_reconstructb,ptsize,cols,'fill','s');xlabel('intercept'),ylabel('gradient');title('hodogram with no wavelet stretch')
% end


% -----------if wavelet stretches across offset------------------
r_synthetic_s=zeros((2*n_s_h),n_fold);

counter=1;

for i=st_vol:inc_vol:end_vol;
    %[wavelet_s,t]=ricker(freq_c_s(i),n_s_h,job_meta.s_rate/1e6,0,phase);
    
    wavelet_matrix=wavelet_digi.all_wavelets_time{i};
    wavelet_s=wavelet_matrix(2:end,level);  
    wavelet_s= wavelet_s';
    wavelet_s=wavelet_s/max(wavelet_s);
    r_synthetic_s(:,counter)=conv(wavelet_s,r(counter));
%     angles2(counter)=angles(i);
    counter=counter+1;
end
[intercept_reconstruct_s,gradient_reconstruct_s] =ig_reconstruct_linefit(r_synthetic_s,angles);


if plotfig==1
    figure(3);
    subplot(2,2,1);scatter(angles,r);xlabel('angle (degrees)');ylabel('reflectivity');title('AVA of reflectivity');
    subplot(2,2,3);imagesc(r_synthetic_s);xlabel('angle (degrees)');ylabel('time sample');title('Synthetic AVA gather')
    subplot(2,2,2);plot(angles,r_synthetic_s');title('Reconstructed AVA curves')
    subplot(2,2,4);scatter(intercept_reconstruct_s,gradient_reconstruct_s);xlabel('intercept'),ylabel('gradient');title('hodogram with real wavelet');
end


% % a plot showing the sample values colouring the hodogram
% if plotfig==1
%     figure(30); 
%     %subplot(2,2,1);scatter(angles,r);xlabel('angle (degrees)');ylabel('reflectivity');title('AVA of reflectivity');
%     subs = 45;
%     ends = 85;
%     r_syntheticb = r_synthetic_s(subs:ends,:);
%     subplot(2,5,[1 6]);imagesc(r_syntheticb);xlabel('angle (degrees)');ylabel('time sample');title('Synthetic AVA gather using extracted wavelets')
%     subs = 45;
%     ends = 65;
%     r_syntheticb = r_synthetic_s(subs:ends,:);
%     intercept_reconstructb = intercept_reconstruct_s(subs:ends);
%     gradient_reconstructb = gradient_reconstruct_s(subs:ends);
%     subplot(2,5,2);imagesc(r_syntheticb);xlabel('angle (degrees)');ylabel('time sample');title('Top half of Synthetic AVA gather')
%     %samples = repmat((1:size(r_synthetic,1))',1,size(r_synthetic,2));
%     
%     samples = repmat((1:size(r_syntheticb,1))',1,10);
%     subplot(2,5,3);imagesc(samples);xlabel('');ylabel('time sample');title('sample colour')
%     ptsize = 60;
%     cols = linspace(1,10,length(intercept_reconstructb));
%     subplot(2,5,4:5);scatter(intercept_reconstructb,gradient_reconstructb,ptsize,cols,'fill','s');xlabel('intercept'),ylabel('gradient');title('hodogram of IG from top half of wavelet')
%     
%     subs = 66;
%     ends = 85;
%     r_syntheticb = r_synthetic_s(subs:ends,:);
%     intercept_reconstructb = intercept_reconstruct_s(subs:ends);
%     gradient_reconstructb = gradient_reconstruct_s(subs:ends);
%     
%     subplot(2,5,7);imagesc(r_syntheticb);xlabel('angle (degrees)');ylabel('time sample');title('Lower half Synthetic AVA gather')
%     %samples = repmat((1:size(r_synthetic,1))',1,size(r_synthetic,2));
%     
%     samples = repmat((1:size(r_syntheticb,1))',1,10);
%     subplot(2,5,8);imagesc(samples);xlabel('');ylabel('time sample');title('sample colour')
%     ptsize = 60;
%     cols = linspace(1,10,length(intercept_reconstructb));
%     subplot(2,5,9:10);scatter(intercept_reconstructb,gradient_reconstructb,ptsize,cols,'fill','s');xlabel('intercept'),ylabel('gradient');title('hodogram of IG from lower half of wavelet')    
%     %colormap(jet);
%     colormap(gray);
% end

% a plot showing the sample values colouring the hodogram
% if plotfig==1
%     figure(40); 
%     %colcj = [0 0 (143/255); (128/255) 1 (128/255); 0 0 (143/255)];
%     %[ cjrgb ] = make_colormap( [0 0.75 2 3 4 5 6 7 8],[0 0 128 255 255 255 128 0 0],[0 0 255  255 0 255 255 0 0],[143 255 128  0 0 0 128 255 143],129,255 );
%     %[ cjrgb ] = make_colormap( [0 1 2  3  4 5 6],[0 0 128 255  128 0 0],[0 0 255   0  255 0 0],[143 255 128   0  128 255 143],129,255 );
%     [ cjrgb ] = make_colormap( [0 1 2  3  4 5 6],[0 0 64 255  64 0 0],[0 0 255   0  255 0 0],[143 255 64   0  64 255 143],127,255 );
%     %colormap(jet)
%     colormap(cjrgb)
%     %subplot(2,2,1);scatter(angles,r);xlabel('angle (degrees)');ylabel('reflectivity');title('AVA of reflectivity');
%     subs = 45;
%     ends = 85;
% 
%     r_syntheticb = r_synthetic_s(subs:ends,:);
%     intercept_reconstructb = intercept_reconstruct_s(subs:ends);
%     gradient_reconstructb = gradient_reconstruct_s(subs:ends);
%     timeb = t(subs:ends);
%     
%     samples = repmat((1:size(r_syntheticb,1))',1,50);
%     
%     
%     %subplot(1,5,2);plot(gradient_reconstructb,timeb);title('Gradient with real wavelet');xlim([-2*abs(gradient),2*abs(gradient)]);grid('on');   
%     xshifts = (150:175);
%     subplot(1,3,1);imagesc(xshifts,timeb,samples);xlabel('');ylabel('time sample');title('sample colour')
%     hold on;
%     plot(intercept_reconstructb,timeb,'-r');title ('Intercept using extracted wavelets','FontSize',14);xlim([-2*abs(intercept),2*abs(intercept)]);grid('on');
%     hold off;
% 
%     xshifts = (300:350);
%     subplot(1,3,2);imagesc(xshifts,timeb,samples);xlabel('');ylabel('time sample');title('sample colour')
%     hold on;
%     plot(gradient_reconstructb,timeb);title('Gradient using extracted wavelets','FontSize',14);xlim([-2*abs(gradient),2*abs(gradient)]);grid('on');
%     %set(p,'Color','black','LineWidth',5)
%     hold off;
%     ptsize = 80;
%     cols = linspace(1,10,length(intercept_reconstructb));
%     subplot(1,3,3);scatter(intercept_reconstructb,gradient_reconstructb,ptsize,cols,'fill','s');xlabel('intercept'),ylabel('gradient');title('IG crossplot showing hodogram with wavelet stretch','FontSize',14);grid('on')
%     xlim([-40 120]); ylim([-100 250]);
% end
% 
% 
% 
% if plotfig==1
%     figure(33);
%     %colormap('bone')
%     subplot(1,5,1);scatter(angles,r);xlabel('angle (degrees)');ylabel('reflectivity');title('synthetic reflectivity');
%     subplot(1,5,2:3);seisplot(r_syntheticb,20,timeb,angles_plot);xlabel('angle (degrees)');ylabel('time sample');title('Synthetic AVA gather using angle variant extracted wavelets');grid('on');
%     %subplot(1,4,2);imagesc(r_syntheticb);xlabel('angle (degrees)');ylabel('time sample');title('Synthetic AVA gather using angle variant extracted wavelets')
%     xlim([0 ((angles_plot(end)*20)+200)]);
%     %set(gca,'XTick',-pi:pi/2:pi)
%     set(gca,'XTickLabel',st_ang:5:end_ang+20);
%     subplot(1,5,4);plot(intercept_reconstructb,timeb,'-r');title ('Intercept from synthetic gather');xlim([-2*abs(intercept),2*abs(intercept)]);grid('on');
%     subplot(1,5,5);plot(gradient_reconstructb,timeb);title('Gradient from synthetic gather');xlim([-2*abs(gradient),2*abs(gradient)]);grid('on');
% end


%--------- plot diferent versions of intercept and gradient traces ----------



if plotfig==1
    figure(4);
    subplot(2,3,1); plot(intercept_orig,t,'-r');title ('Original Intercept');xlim([-2*abs(intercept),2*abs(intercept)]);grid('on');
    subplot(2,3,2);plot(spline(t,intercept_reconstruct,t_fine),t_fine,'-r');title ('Intercept without stretch');xlim([-2*abs(intercept),2*abs(intercept)]);grid('on');
    subplot(2,3,5);plot(spline(t,gradient_reconstruct,t_fine),t_fine);title ('Gradient without stretch');xlim([-2*abs(gradient),2*abs(gradient)]);grid('on');
    subplot(2,3,4); plot(gradient_orig,t);title('Original Gradient');xlim([-2*abs(gradient),2*abs(gradient)]);grid('on');
    subplot(2,3,3);plot(spline(t,intercept_reconstruct_s,t_fine),t_fine,'-r');title ('Intercept with real wavlet');xlim([-2*abs(intercept),2*abs(intercept)]);grid('on');
    subplot(2,3,6);plot(spline(t,gradient_reconstruct_s,t_fine),t_fine);title('Gradient with real wavelet');xlim([-2*abs(gradient),2*abs(gradient)]);grid('on');
end

% figure
% for ii = 2:size(r_syntheticb(:,2:5),2)
%     seisplot(r_syntheticb(:,ii),150)
% end

end

function seisplot(x1,x2,x3,ang)
    n=size(x1,1); z1=zeros(1,n);
% x1=x1/max(abs(x1)); % Optional
    %yt=0; % Base line = 0
    yt = 0;
    for ii = 1:size(x1,2)    
        z1(1:n)=x1(:,ii)+yt;
        for i=1:n, w1(i)=yt; if z1(i) > yt, w1(i)=z1(i); end; end;
        w1(1)=yt; w1(n)=yt; x=1:n;
        %plot(x,w1,'-k'); patch(x,w1,'k'); 
        hold; plot(z1,x3,'-k'); hold off;
        yt = yt + x2;
    end
end

function [s,t] = ricker(f,n,dt,t0,phase)
%RICKER creates an causal ricker wavelet signal

ar =(-n+0.5):1:(n-0.5);
t = ar*dt;
tau = t-t0;
%s = (1-tau.*tau*f^2*pi^2).*exp(-tau.^2*pi^2*f^2+phase);
s = (1-tau.*tau*f^2*pi^2).*exp(-(tau*pi*f-phase).^2);
end

%function to reconstruct intercept and gradient by line fitting
function [int,grad]=ig_reconstruct_linefit(synthetic,angles)

angles = angles(~isnan(mean(synthetic)));
synthetic = synthetic (:,~isnan(mean(synthetic)));
angle_t=(sind(angles).*sind(angles));
int=zeros(1,size(synthetic,1));
grad=zeros(1,size(synthetic,1));
for n_t=1:size(synthetic,1)
    reflections=synthetic(n_t,:);
    p=polyfit(angle_t,reflections,1);
    int(n_t)=p(2);
    grad(n_t)=p(1);
end

end
function [ rgb ] = make_colormap( color_idx,red,green,blue,num_colors,max_color )
% make a colormap from rgb values
%   color_idx: index of supplied rgb values (vector)
%   red: red values 
%   green: green values 
%   blue: blue values 
%   num_colors: number of colours in output colormap
%   max_color: [max_color max_color max_color] is white
%   
%   color_idx, red, green and blue vectors should all be same size
%   color_idx is zero-indexed
%   
%   Example:
%   red_white_blue = make_colormap([0 1 2],[1 1 0],[0 1 0],[0 1 1],256,1);

color_idx = color_idx./(max(color_idx)/num_colors);

num_colors = num_colors-1;

red_interp = interp1(color_idx,red,[0:num_colors])';
green_interp = interp1(color_idx,green,[0:num_colors])';
blue_interp = interp1(color_idx,blue,[0:num_colors])';

rgb = [red_interp green_interp blue_interp];
rgb = rgb./max_color;

end
