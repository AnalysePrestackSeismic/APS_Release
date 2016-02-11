function [output] =ava_wedge_digi(vp1,vs1,rho1,vp2,vs2,rho2,job_meta_path,wavelet_file,z,method,plot_results)
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

%% Function to investigae eva modeling by modeling a wedge rather than a single interface. 
% This  helps in evaluatin the effect of rock properties along with the
% effect of tuning on the Ampltiude variation with Offset.
% INPUTS:
%   vp1 = p-wave velocity of first interface (background),
%   vs1 = s-wave velocity of first interface (background),
%   rho1 = density of first interface (background),
%   vp2 = p-wave velocity of reservoir body ,
%   vs2 = s-wave velocity of reservoir body ,
%   rho2 = density velocity of reservoir body,
%   max_angle = Maximum Angle in degrees you want to model,
%   freq_c = Central Frequency of wavelet,
%   method = 1 : Aki Richar'd 3 Term, 2: Zeopritz, 3: 2 Term, Shuey
%   plot_results=Make this 1 to plot stuff


                                                 % Cloas all open figures
job_meta = load(job_meta_path);
% Create the wavelet, currently no stretching is introduced

wavelet_path = strcat(job_meta.wav_directory,wavelet_file);
wavelet_digi=load(wavelet_path);

st_ang=floor(job_meta.tkey_min);
inc_ang=floor(job_meta.tkey_inc);
end_ang =floor(job_meta.tkey_max);
st_vol=1;
inc_vol=1;
end_vol=1+(end_ang-st_ang)/inc_ang;
t_wavelet=(0:127)*job_meta.s_rate/1e6;


%% Create the wedge model / set options for geometries 
wedge_model_thickness = 0.1;                               % Define Wedge Model max thickness in s
v_sampling=0.00025;%job_meta.s_rate/1e6;                                         % Define Wedge Modeling vertical Sampling in s (Keep this fairly low value so that figures are not blocky)
h_sampling=1;                                               % Horizontal Sampling defined as 1 m. This doent matter much to the actual results (since its in thickness
v_padding=0.5*wedge_model_thickness/v_sampling;             % Vertical padding
n_s_h=1/v_sampling;                                         % Number of traces (or horizonatl samples)
n_v=2*v_padding+1+wedge_model_thickness/v_sampling;         % Total number of samples  vertically including padding
n_h=1+wedge_model_thickness/v_sampling;                     % Total number of samples horizontally
n_fold= end_vol;                                           % Number of angles

% Define axes
t_axis =(1:n_v)* v_sampling*1000;                           % Create time axis in ms
h_axis =(1:n_h)*h_sampling;                                 % Create horizontal axis in ms

wedge_model_mat=zeros(n_v,n_h);                             % Pre allocate memory to Wedge Model
wedge_model_mat(v_padding,:) = 1;                           % Create Top of the Wedge Model as +1
th=v_padding+1;                                             % Initialize minimum thickness as 1 sample
wedge_model_r=zeros(n_v,n_h,n_fold);                           % Pre allocate memory for storing AVA reflectivity gather
wedge_model_synth=zeros(n_v,n_h,n_fold);                       % Pre allocate memory for storing AVA synthetic gather (post convolution )                     

% Create the Base of the wedge model ., Slanted base
for i=1:n_h
    wedge_model_mat(th,i)=-1;                               % Create Top of the Wedge Base as +1               
    th=th+1;                                                % Increment thickness in samples . 1 sample per loop iteration
end

body=cumsum(wedge_model_mat);                               % Make this into a body like thing by intergration
time_thickness_true=(1:n_h)*v_sampling*1000;                % Vector of true time thickness of the wedgmodel for every trace
%% Modeling




if method ==1
    [~,r_ava]=avo_3term_aki_richard(vp1,vs1,rho1,vp2,vs2,rho2,end_ang); % Create reflectivity using ava modeling equation
end
if method ==2
    [~,r_ava]=avo_zeoppritz(vp1,vs1,rho1,vp2,vs2,rho2,end_ang); % Create reflectivity using ava modeling equation
end

ang=st_ang:inc_ang:end_ang;
r_ava=r_ava(ang);

% Create 2D sythtetic elatic reflectivity volumes ( Angle being 3rd
% diamension anfd then produce stnthetic gather
% tuning_curve(wavelet,1);
% -----------if wavelet stretches across offset------------------
% r_synthetic_s=wedge_model_r*0;


counter=0;
level=floor(z/(32*job_meta.s_rate/1e3));                                            % index level from which wavelet will be taken
for i=st_vol:inc_vol:end_vol;
    wedge_model_r(:,:,i)=wedge_model_mat*r_ava((i));
    
    wavelet_matrix=wavelet_digi.all_wavelets_time{i};
    wavelet_s=wavelet_matrix(2:end,level);
    wavelet_s=wavelet_s/max(wavelet_s);
    
    % Initialize matrix for storing wavelets
    if i==st_vol
        wavelet_set=zeros(length(wavelet_s),n_fold);
    end
    
    % Check if wavlet is empty or not
    if sum(isnan(wavelet_s))==0
        counter=counter+1;
        wavelet_set(:,i)=wavelet_s;
        wavelet_interp=interp1(t_wavelet,wavelet_s,0:v_sampling:t_wavelet(end),'spline');
        wedge_model_synth(:,:,counter)=conv2(squeeze(wedge_model_r(:,:,i)),wavelet_interp','same');
        
    else
        ang(i)=-1;
    end    
    
end
ang=ang(ang>-1);
r_ava=r_ava(ang>-1);
wedge_model_synth=wedge_model_synth(:,:,1:counter);


%% Calculate 2D AVA

%   Pick the top and base rflections on the gather
ava_refl_top=squeeze(min( wedge_model_synth,[],1));
ava_refl_base=squeeze(max( wedge_model_synth,[],1));
int=zeros(1,n_h);
grad=zeros(1,n_h);
chi_angle=zeros(1,n_h);
       
intercept=wedge_model_synth(:,:,1);                         % Intercept is the nearest traces
clear wedge_model_synth wedge_model_r;

% Calculate intercept and Gradient and the Top horizon by Linear regression
% through the gather

for ii=1:n_h
    p=polyfit((sin(ang(ang<35)*pi/180)).^2,ava_refl_top(ii,(ang<35)),1);   
    int (ii)=p(2);
    grad(ii)=p(1);
    chi_angle(ii)=chi_transform(int (ii),grad(ii));
end

 
% Pick the apparent tops and base on the nearest stack or intercept
% intercept = wedge_model_synth(:,:,1);
[horizon_top_amp,horizon_top_idx] =max(intercept);
[horizon_base_amp,horizon_base_idx]=min(intercept);
time_thickness_ap = abs((horizon_base_idx -horizon_top_idx))*v_sampling*1000;
bli=cumsum(intercept); % Create band limited impedance by trace integration

output =struct;
output.ang=ang;
output.r_ava=r_ava;
output.ava_refl_top=ava_refl_top;
output.time_thickness_true=time_thickness_true;
output.time_thickness_ap=time_thickness_ap;
output.int=int;
output.grad=grad;
output.t_axis=t_axis;
output.h_axis=h_axis;
output.t_wavelet=t_wavelet;
output.wavelet_set=wavelet_set;
output.chi_angle=chi_angle;



%% Plot Results
if plot_results ==1
    close all; 
    red_white_blue = make_colormap([0 1 2],[1 1 0],[0 1 0],[0 1 1],1024,1);
    %---------------- FIGURE : 1 ------------------------------------------
    
    % Plot the Wedge Model
    figure(1);
    subplot(1,4,1); 
    imagesc(body);colormap('jet');
    title('Wedge reflectivity model');ylabel('Time (ms)');xlabel('Horizontal Dist (m)')
    set(gca,'YTick',floor(n_v/10):floor(n_v/10):n_v);
    set(gca,'YTickLabel',t_axis(floor(n_v/10):floor(n_v/10):n_v));
    
    % Plot Synthetic Seismogram from Wedge Model
    subplot(1,4,2); 
    imagesc(intercept);colormap(red_white_blue);
    title('Synthetic Intercept');ylabel('Time (ms)');xlabel('Horizontal Dist (m)')
    hold on; 
    plot(horizon_top_idx,'--g','LineWidth',3);
    plot(horizon_base_idx,'--k','LineWidth',3);hold off;
    set(gca,'YTick',floor(n_v/10):floor(n_v/10):n_v);
    set(gca,'YTickLabel',t_axis(floor(n_v/10):floor(n_v/10):n_v));
    
    % Plot Band Limited Impedance
    subplot(1,4,3); 
    imagesc(bli);
    title('Band Limited Impedance');colormap(red_white_blue);
    ylabel('Time (ms)');xlabel('Horizontal Dist (m)')
    hold on; 
    plot(horizon_top_idx,'--g','LineWidth',3);
    plot(horizon_base_idx,'--k','LineWidth',3);hold off;
    set(gca,'YTick',floor(n_v/10):floor(n_v/10):n_v);
    set(gca,'YTickLabel',t_axis(floor(n_v/10):floor(n_v/10):n_v));
       
    % Plot tuning curves( Amplitude vs Aparent Thickness) for the top and     % base
    subplot(1,4,4);
    plot([0 max(time_thickness_ap)],[0 0],'g','LineWidth',2);hold on;
    plot(time_thickness_ap,horizon_top_amp,'-r','LineWidth',3);hold on;
    plot(time_thickness_ap,horizon_base_amp,'-b','LineWidth',3);hold off;
    xlim([0  max(time_thickness_ap)]);
    xlabel('Apparent Time Thickness (ms)');ylabel('Amplitude'); title('Tuning Curve');
    
    %------------------- FIGURE: 2 --------------------------------------------------
    % Plot Reflectivity vs angle from half space model
    figure(2)
    subplot(3,2,1);
    plot(ang,r_ava,'-g','LineWidth',3);hold on;
   plot([0 max(ang)],[0 0],'k','LineWidth',2);hold off;
    xlabel('Angle');ylabel('Reflectivity');title('ava plot');
    ylim([-1 1]);
    
    % Plot wavelet
    subplot(3,2,2);
    plot(t_wavelet,wavelet_set,'LineWidth',2);
    xlabel('time');title('Wavlelet');
    xlim([0 max(t_wavelet)]);
%     set(gca,'YTick',20:120);
%     set(gca,'YTickLabel',t_wavelet(20:120));
    
    % For  Top: Plot Reflectivity (as colour) with changing angle and wedge thickness
    subplot(3,2,3);
    imagesc(ava_refl_top);
    n=size(ava_refl_top,1);
    set(gca,'YTick',floor(n/10):floor(n/10):n);
    set(gca,'YTickLabel',time_thickness_true(floor(n/10):floor(n/10):n));
    set(gca,'XTick',1:size(ava_refl_top,2));
    set(gca,'XTickLabel',ang(1:size(ava_refl_top,2)));
    
    ylabel('True Time Thickness (ms)');   
    title('AVA Top');xlabel('Angle');
    
%     % For Base: Plot Reflectivity (as colour) with changing angle and wedge thickness
%     subplot(3,2,4);
%     imagesc(ava_refl_base);
%     n=size(ava_refl_base,1);
%     set(gca,'YTick',floor(n/10):floor(n/10):n);
%     set(gca,'YTickLabel',time_thickness_true(floor(n/10):floor(n/10):n));
%     ylabel('True Time Thickness (ms)');
%     title('AVA Base');xlabel('Angle');
    
    % Plot Intercept and Gradient bs Apparent Time Thickness
    subplot(3,2,4); 
    plot(ang,ava_refl_top);hold on;
   
    plot([0 max(ang)],[0 0],'k','LineWidth',2);
     ylim([-1 1]);hold off;
%     set(gca,'XTick',1:size(ava_refl_top,2));
%     set(gca,'XTickLabel',ang(1:size(ava_refl_top,2)));
    
%     plot(time_thickness_ap,int,'LineWidth',3);hold on;
%     plot(time_thickness_ap,grad,'-r','LineWidth',3);
%     plot([0 max(time_thickness_ap)],[0 0],'g','LineWidth',2);hold on;
%     xlim([0 max(abs(time_thickness_true))]);hold off; 
%     title('Intercept(blue), Gradient(red)');xlabel('Aparent Time Thickness (ms)'); ylabel('Amplitude')
    
    % Plot Intercept Gradient Crossplot for the Wedge.
    subplot(3,2,5); 
    scatter(int,grad,5,[0 1 1]);
    xlim([-max(abs(int)) max(abs(int))]);
    ylim([-max(abs(grad)) max(abs(grad))]);
    grid on;
    xlabel('intercept');ylabel('Gradient');title('Tuning IG Hodogram');
    
    % Plot chi angle with true thickness
    subplot(3,2,6); 
    scatter(time_thickness_true,chi_angle,5,[0 0 1]);
    xlim([0 max(abs(time_thickness_true))]);
    ylim([0 360]); xlabel('True time thickness');
    ylabel('Polar_Angle(I-G)');
    
    %-------------------------------------------------------------------- 
end

end
function [ rgb ] = make_colormap( color_idx,red,green,blue,num_colors,max_color )
%%  make a colormap from rgb values
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
function[chi_ang]=chi_transform(intercept,gradient)
%function to transform an intercept gradient pair into chi angle
quad=2*(intercept>0)+(gradient>0);
ang1=atand(gradient/intercept);

%Unwrap the angles
switch quad
    case 0
        chi_ang =270+ang1;
    case 1
        chi_ang = 270+ang1;
    case 2
        chi_ang = 90+ang1;
    case 4
        chi_ang = 90+ang1;
end

end
