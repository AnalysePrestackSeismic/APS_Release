function [ava_refl_top,time_thickness_true,time_thickness_ap] =ava_wedge(vp1,vs1,rho1,vp2,vs2,rho2,max_angle,freq_c,method,plot_results)

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


close all;                                                  % Cloas all open figures

%% Create the wedge model / set options for geometries 
wedge_model_thickness = 0.06;                               % Define Wedge Model max thickness in s
v_sampling=0.00005;                                         % Define Wedge Modeling vertical Sampling in s (Keep this fairly low value so that figures are not blocky)
h_sampling=1;                                               % Horizontal Sampling defined as 1 m. This doent matter much to the actual results (since its in thickness
v_padding=0.5*wedge_model_thickness/v_sampling;             % Vertical padding
n_s_h=1/v_sampling;                                         % Number of traces (or horizonatl samples)
n_v=2*v_padding+1+wedge_model_thickness/v_sampling;         % Total number of samples  vertically including padding
n_h=1+wedge_model_thickness/v_sampling;                     % Total number of samples horizontally
n_a= 1+max_angle;                                           % Number of angles

% Define axes
t_axix =(1:n_v)* v_sampling*1000;                           % Create time axis in ms
h_axis =(1:n_h)*h_sampling;                                 % Create horizontal axis in ms


wedge_model_mat=zeros(n_v,n_h);                             % Pre allocate memory to Wedge Model
wedge_model_mat(v_padding,:) = 1;                           % Create Top of the Wedge Model as +1
th=v_padding+1;                                             % Initialize minimum thickness as 1 sample
wedge_model_r=zeros(n_v,n_h,n_a);                           % Pre allocate memory for storing AVA reflectivity gather
wedge_model_synth=zeros(n_v,n_h,n_a);                       % Pre allocate memory for storing AVA synthetic gather (post convolution )                     

% Create the Base of the wedge model ., Slanted base
for i=1:n_h
    wedge_model_mat(th,i)=-1;                               % Create Top of the Wedge Base as +1               
    th=th+1;                                                % Increment thickness in samples . 1 sample per loop iteration
end

body=cumsum(wedge_model_mat);                               % Make this into a body like thing by intergration
time_thickness_true=(1:n_h)*v_sampling*1000;                % Vector of true time thickness of the wedgmodel for every trace
%% Modeling

% Create the wavelet, currently no stretching is introduced
[wavelet,t]=ricker(freq_c,n_s_h,v_sampling,0,0);            % create wavelet (only ricker right now)
wavelet=-1*wavelet;                                         % Reverse polarity of wavelet
if method ==1
[ang,r_ava]=avo_3term_aki_richard(vp1,vs1,rho1,vp2,vs2,rho2,max_angle); % Create reflectivity using ava modeling equation
end
if method ==2
[ang,r_ava]=avo_zeoppritz(vp1,vs1,rho1,vp2,vs2,rho2,max_angle); % Create reflectivity using ava modeling equation
end
% r_min=  abs(min(r_ava));
% r_max=  abs(max(r_ava));
% Create 2D sythtetic elatic reflectivity volumes ( Angle being 3rd
% diamension anfd then produce stnthetic gather
% tuning_curve(wavelet,1);

for ii=1:(max_angle+1)
    wedge_model_r(:,:,ii)=wedge_model_mat*r_ava((ii));
    
    wedge_model_synth(:,:,ii)=conv2(wedge_model_r(:,:,ii),wavelet','same');
    
end




%% Calculate 2D AVA

%   Pick the top and base rflections on the gather
ava_refl_top=squeeze(min( wedge_model_synth,[],1));
ava_refl_base=squeeze(max( wedge_model_synth,[],1));
int=zeros(1,n_h);grad=zeros(1,n_h);

% gradient=zeros(n_v,n_h);
% intercept=zeros(n_v,n_h);
% for j=1:n_v
%     for k=1:n_h
%         p=polyfit((sin(ang*pi/180)).^2,[squeeze(wedge_model_synth(j,k,:))]',1);
%          intercept(j,k)=p(2);
%         gradient(j,k)=p(1);
%     end
% end
        
intercept=wedge_model_synth(:,:,1);                         % Intercept is the nearest traces

% Calculate intercept and Gradient and the Top horizon by Linear regression
% through the gather

for ii=1:n_h
    p=polyfit((sin(ang*pi/180)).^2,ava_refl_top(ii,:),1);   
    int (ii)=p(2);
    grad(ii)=p(1);
end
 
% Pick the apparent tops and base on the nearest stack or intercept
% intercept = wedge_model_synth(:,:,1);
[horizon_top_amp,horizon_top_idx] =max(intercept);
[horizon_base_amp,horizon_base_idx]=min(intercept);
time_thickness_ap = abs((horizon_base_idx -horizon_top_idx))*v_sampling*1000;
bli=cumsum(intercept); % Create band limited impedance by trace integration

%% Plot Results
if plot_results ==1
    red_white_blue = make_colormap([0 1 2],[1 1 0],[0 1 0],[0 1 1],1024,1);
    %---------------- FIGURE : 1 ------------------------------------------
    
    % Plot the Wedge Model
    figure(1);
    subplot(1,4,1); 
    imagesc(body);colormap('jet');
    title('Wedge reflectivity model');ylabel('Time (ms)');xlabel('Horizontal Dist (m)')
    set(gca,'YTick',[500 1000 1500 2000]);
    set(gca,'YTickLabel',[t_axix(500) t_axix(1000) t_axix(1500) t_axix(2000)]);
    
    % Plot Synthetic Seismogram from Wedge Model
    subplot(1,4,2); 
    imagesc(intercept);colormap(red_white_blue);
    title('Synthetic Intercept');ylabel('Time (ms)');xlabel('Horizontal Dist (m)')
    hold on; 
    plot(horizon_top_idx,'--g','LineWidth',3);
    plot(horizon_base_idx,'--k','LineWidth',3);hold off;
    set(gca,'YTick',[500 1000 1500 2000]);
    set(gca,'YTickLabel',[t_axix(500) t_axix(1000) t_axix(1500) t_axix(2000)]);
    % subplot(1,4,4); imagesc(gradient);colormap('jet'); title('Gradient');
    % subplot(1,6,5); imagesc(squeeze(wedge_model_synth(:,100,:)));colormap('jet');title('Gather thin end');
    % subplot(1,6,6); imagesc(squeeze(wedge_model_synth(:,400,:)));colormap('jet');title('Gather thick end');
    
    % Plot Band Limited Impedance
    subplot(1,4,3); 
    imagesc(bli);
    title('Band Limited Impedance');colormap(red_white_blue);
    ylabel('Time (ms)');xlabel('Horizontal Dist (m)')
    hold on; 
    plot(horizon_top_idx,'--g','LineWidth',3);
    plot(horizon_base_idx,'--k','LineWidth',3);hold off;
    set(gca,'YTick',[500 1000 1500 2000]);
    set(gca,'YTickLabel',[t_axix(500) t_axix(1000) t_axix(1500) t_axix(2000)]);
       
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
    plot(ang,r_ava,'-g','LineWidth',3);
    xlabel('Angle');ylabel('Reflectivity');title('ava plot');
    
    % Plot wavelet
    subplot(3,2,2);
    plot(t*1000,wavelet,'-k','LineWidth',3);
    xlabel('time');title('Wavlelet');
    
    % For  Top: Plot Reflectivity (as colour) with changing angle and wedge thickness
    subplot(3,2,3);
    imagesc(ava_refl_top);
    set(gca,'YTick',200:200:1200);
    set(gca,'YTickLabel',time_thickness_true(200:200:1200));
    ylabel('True Time Thickness (ms)');   
    title('AVA Top');xlabel('Angle');
    
    % For Base: Plot Reflectivity (as colour) with changing angle and wedge thickness
    subplot(3,2,4);
    imagesc(ava_refl_base);
    set(gca,'YTick',200:200:1200);
    set(gca,'YTickLabel',time_thickness_true(200:200:1200));
    ylabel('True Time Thickness (ms)');
    title('AVA Base');xlabel('Angle');
    
    % Plot Intercept and Gradient bs Apparent Time Thickness
    subplot(3,2,5); 
    plot(time_thickness_ap,int,'LineWidth',3);hold on;
    xlim ([0  max(time_thickness_ap)])
    plot(time_thickness_ap,grad,'-r','LineWidth',3);hold off; 
    title('Intercept(blue), Gradient(red)');xlabel('Aparent Time Thickness (ms)'); ylabel('Amplitude')
    
    % Plot Intercept Gradient Crossplot for the Wedge.
    subplot(3,2,6); 
    plot(int,grad,'-k','LineWidth',3);
    xlim([-max(abs(int)) max(abs(int))]);
    ylim([-max(abs(grad)) max(abs(grad))]);
    grid on;
    xlabel('intercept');ylabel('Gradient');title('Tuning IG Hodogram');
    %-------------------------------------------------------------------- 
end

end

function [s,t] = ricker(f,n,dt,t0,phase)
%% RICKER creates an causal ricker wavelet signal
T = dt*(n-1);
t = (-1*T):dt:T;
tau = t-t0;
%s = (1-tau.*tau*f^2*pi^2).*exp(-tau.^2*pi^2*f^2+phase);
s = (1-tau.*tau*f^2*pi^2).*exp(-(tau*pi*f-phase).^2);
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
