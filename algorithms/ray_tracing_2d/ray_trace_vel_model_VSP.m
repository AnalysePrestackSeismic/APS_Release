function []=ray_trace_vel_model_VSP(volume_path,il_byte_loc,xl_byte_loc,x_interval,direction)
%   Inputs:
%         volume_path = Path of seismic Volume
%         il_byte_loc = Inline byte location
%         xl_byte_loc = xlbyte location
%         shot_loc = shote location in [il xl]
%         x_interval = inline or xline interval in m (depending on direction)
%         direction = 1 if its an inline  0 if its a xline
%         output_dir = output directory
%         smo = smoothening parameter
%% Parameter
s=2.5;                    % Parameter to control the maximum time propagation (in sec)of the rays
s_a1 = 1 ;             % to control the spread of the ray fan. A bigger number means bigger spread
s_a2= 1000;              % to control the decimation of the ray fan. 
nla=0;                  % Make this 1 for nonlinear spread of angles (density of rays is more when they are more vertical) , anyother value for linear spread    
r_vel = 1545;           % Mimimum rock velocity
% VSP Parameters
%threshold_dist=10;      % Threshold vertical distance in m
rec_spac =15;        % VSP receiver Spacing in meters
option_intersection=1; % Calculate ray's intersection with linear solver if 1 and by averaging if anything else
% Plot parameters
%offset = 9000;         % Maximum Offset in Plots
z_max=6000;             % Maximum Z in Ray Tracing Plot
v_max=3500;             % Maximum Velocity  (to limit color) in Ray Tracing Plot

%% Define the mat file to be used by the ray tracer

close all;                                                                      % Close all open figures
[seis_meta il_xl_pos seismic]=segy_to_mat(il_byte_loc,xl_byte_loc,volume_path); % Function in extras that converts a segy file to a mat file
z = 0:(size(seismic,1)-1);
vel =seismic;                                                                   % The velocity matrix (2D for now)
if direction==1
    x= double(il_xl_pos (:,2)-il_xl_pos (1,2))*x_interval;                      % Define x axis if the line is an inline
elseif direction==0
    x= double(il_xl_pos (:,1)-il_xl_pos (1,1))*x_interval;                      % Define x axis if the line is a xline
else
    fprintf ('\nrespecify right dir: 1 if its an inline  0 if its a xline\n')
end
x=x';                                                                           % Transpose x axis to make a row vector
z=z'*seis_meta.s_rate/1000;                                                     % Define z axis as a column vector

%% Pick water bottom
v_temp = (vel-r_vel);
v_temp(v_temp<0)=9999;
[~,i] = min(v_temp(5:size(v_temp,1),:));
z_wb=z(i);                                                                      % Water bottom depth
clear v_temp;

%% Well Path

well_z_points=[0 1000 2000 3000 4500 5000];                                    % Deviation survey: Well path discrete z values
well_x_points=[10000 10000 10200 11000 11500 13200];                          % Deviation survey: Well path discrete xvalues

well_top_x=min(well_x_points);                                                      % Well top x location
well_top_z=min(well_z_points);                                                      % Well top z location

well_bot_x=max(well_x_points);                                                      % Well bottom x location
well_bot_z=max(well_z_points);                                                      % Well bottom z location

well_path_zz=z(z>well_top_z & z<well_bot_z);                                        % Crop the well path z values, grid of seismic velcoity is used
well_path_xx=interp1(well_z_points,well_x_points,well_path_zz,'linear','extrap');   % Resample the well path x location to the corresponding z values
well_path_vel=zeros(size(well_path_zz));

% Extract the velocity along the well path
for p=1:length(well_path_vel)
    well_path_vel(p)=vel(well_path_zz(p)/(seis_meta.s_rate/1000),round((well_path_xx(p)-x(1))/x_interval)+1);
end

%% -------------------------Install Velocity Model--------------------------------------------------------
nx=length(x);
nz=length(z);
dx=x(2)-x(1);
dg=z(2)-z(1);

rayvelmod_bg(vel,dx,dg);                                            % Install the velocity model
close all;                                                          % Close all open figures
figure(1);
imagesc(x,z,vel,[1400 v_max]);ylim([0 z_max]);                      % Plot the velocity model
hold all;
%plot(x,z_wb,'-k','LineWidth',2);                                   % Plot water bottom
plot(well_path_xx,well_path_zz,'-k','LineWidth',3);                 % Plot well path
%% ------------------Estimate tmax,dt,tstep----------------------
vlow=min(min(vel));
tmax=s*max(z)/vlow;dt=.004;tstep=0:dt:tmax;                         % Define the time step vector
%% ----------------Define Ray Fan-----------------------------------

angles=(-pi/2+pi/(2.1+s_a1)):pi/s_a2:(pi/2-pi/(2.1+s_a1));                      % Define the angles controling the spread of ray fan
% For non-linearly decimated ray spread
if nla==1
    sclr_l=abs(-1:2/(length(angles)-1):1);                          % Linear scalar vector between 1 and 0
    sclr_nl=(sclr_l).^(1.7);                                        % make  vector non linear using a power, increase power for more non linearity
    angles=((-pi/2+pi/s_a1):pi/s_a2:(pi/2-pi/s_a1)).*sclr_nl;       % Scale Angle vector non linearly. Max and Min angle remmains constant.
end
%% ----------get input of shot location------------------------------------------
x0=round(nx/3)*dx;z0=0;                                             % Define the default location of shot

disp(['Source currently set at x=',num2str(x0),' z0=',num2str(z0)])
msgx='Enter x coordinate of new source location (<cr> for no change : )';
xxxx=input(msgx);
if(~isempty(xxxx))
    x0=xxxx;
end
msgz='Enter z coordinate of new source location (<cr> for no change :)';
zzzz=input(msgz);
if(~isempty(zzzz))
    z0=zzzz;
end

%% -----------------Ray Tracing------------------------------------------------
indx=near(x,x0);indz=near(z,z0);
v0=vel(indz,indx);
sr=zeros(length(angles),6);                             % Matrix to keep shot record
td_pairs=zeros(length(angles)*2,3);                     % Matrix to keep t-d pairs

t1=clock;
k_min=20;
counter=0;
for k=1:length(angles)
    r0=[x0 z0 sin(angles(k))/v0 cos(angles(k))/v0];     % Define the initial ray vector
    [t,r]=shootrayvxz_bg(tstep,r0);                     % Shoot the rays
    %line(r(:,1),r(:,2),ones(size(t)),'color','w');     % Plot the rays
    sr(k,1)=angles(k);                                  % Store angles a first column
    z_k=r(:,2);                                         % This is the array of z coordinates in the ray's travel path
    x_k=r(:,1);                                         % This is the array of x-coordintes in the ray's travel path
    [z_min ind1]=min(abs(z_k(k_min:length(z_k))));      % Remove first couple of points near datum and Find when the ray reaches the datum
    [sr(k,4),ind_maxz1]=max(z_k(1:length(z_k)));        % Store the the maximum z the ray sees
    
    % If the ray reaches the datum at all
    if (z_min<z_k(k_min-1))
        sr(k,2)=r((k_min+ind1-1),1);                    % Store the x location at which the ray is received
        sr(k,3)=t(k_min+ind1-1);                        % Store the the total travel time when the ray is received
        
        [c,wb_x_index]=min(x-r(ind_maxz1,1));           % Find index of x-axis at the location of maximum depth reached
        wb_z_k=z_wb(wb_x_index);                        % Water bottom depth correspoinding to maximum depth reached
        sr(k,6)=sr(k,4);
        sr(k,5)=sr(k,4)-wb_z_k;                         % the overburden thickness corresponding to the maximum depth reached
    end
    
    %----------Checkshot Creation----------------
    [~,ind_maxz2]=min(abs(well_path_zz-sr(k,4)));                      % Find the deepest z location in the ray path in the original zgrid
    
    wellz_resample=interp1(well_path_xx(1:ind_maxz2),well_path_zz(1:ind_maxz2),x_k(1:ind_maxz1),'linear','extrap'); % Resample the well path on the ray_trace grid
    
    del_z=z_k(1:ind_maxz1)-wellz_resample;
   
    
    % Find intersections, even multiple ones
    for j=2:ind_maxz1
        % check if any sign change occurs and this doent fall in
        % extrapolated section of well
        if(del_z(j)*del_z(j-1)<0 && wellz_resample(j)<well_bot_z)
            counter=counter+1;                                           % Increment Counter, first location in t-d pair remains initialized at t=0 d=0
            if option_intersection==1
                [td_pairs(counter,1) td_pairs(counter,3)] =linear_eq_solve(z_k(j),x_k(j),z_k(j-1),x_k(j-1),wellz_resample(j),x_k(j),wellz_resample(j-1),x_k(j-1)); % Calcutae z value of TD pair
                
                td_pairs(counter,2)=t(j-1)+((t(j)-t(j-1))*(td_pairs(counter,1)-z_k(j-1))/(z_k(j)-z_k(j-1)));                         % Calcutae z value of TD pair
                x_k(j)=td_pairs(counter,3);
                z_k(j)=td_pairs(counter,1);                
            else
                
                td_pairs(counter,1)=(z_k(j)+z_k(j-1))/2;                     % Calcutae z value of TD pair
                td_pairs(counter,2)=(t(j)+t(j-1))/2;                         % Calcutae z value of TD pair
                td_pairs(counter,3)=(x_k(j)+x_k(j-1))/2;                         % Calcutae x value of TD pair
                x_k(j)=(x_k(j)+x_k(j-1))/2;
                z_k(j)=(z_k(j)+z_k(j-1))/2;                
            end
            line(x_k(1:j),z_k(1:j),ones(j,1),'color','w');               % Plot the rays uptil they reach well
        end
    end
end


td_pairs=td_pairs(1:counter,:);                                                 % Remove the unaltered zeros values of the pre-initialized t-d pairs matrix
clear c c1 j intersection wb_z_k wb_x_index z_min ind1 ind_maxz2 ind_maxz1 x_k z_k distz_well ind2 wellz_resample k i r r0 r_vel t counter wellz_resample_crop z_k_crop z_k_crop2 t_crop t_crop2;
title('Ray Tracing','FontSize',25); xlabel('x axis (m)','FontSize',15); ylabel('z axis (m)','FontSize',15);
set(gca,'XAxisLocation','top');ylim([0 z_max]);


well_z_cs = well_top_z:rec_spac:well_bot_z;                                     % Check Shot Aquisition Z axis
well_t_cs = interp1(td_pairs(:,1),td_pairs(:,2),well_z_cs,'spline','extrap');   % Check Shot Aquisition Time Values resampled
well_t_cs(well_t_cs<min(td_pairs(:,2)))=min(td_pairs(:,2));

well_x_cs = interp1(td_pairs(:,1),td_pairs(:,3),well_z_cs,'spline','extrap');


   
%Raw  velcoity from depth and first arrivals
vel_cs_raw = diff(well_z_cs )./diff(well_t_cs);
vel_cs_raw = [ vel_cs_raw(1) vel_cs_raw];

%Corrected Velocity
well_t_cs_corr = well_t_cs.*(well_z_cs./((well_z_cs.^2+(well_x_cs-x0).^2).^0.5)); % Correction of travel time for offset
vel_cs_corr = diff(well_z_cs )./diff(well_t_cs_corr);
vel_cs_corr = [ vel_cs_corr(1) vel_cs_corr];

%% ----------------Plot VSP-------------------------------------
figure(2);
subplot(1,3,1);
scatter(td_pairs(:,1),td_pairs(:,2)*1000,15,'filled','MarkerEdgeColor','b','MarkerFaceColor','c');ylim([0 max(td_pairs(:,2))*1000*1.15]);
title('Ray Traced Check shot','FontSize',25);xlabel('Depth(m)','FontSize',15);
ylabel('Arrival time(ms)','FontSize',15);
grid(gca,'minor');
set(gca,'YDir','reverse');


subplot(1,3,2);
scatter(well_z_cs,well_t_cs*1000,15,'filled','MarkerEdgeColor','b','MarkerFaceColor','c');ylim([0 max(td_pairs(:,2))*1000*1.15]);
title('Check shot(at aquired res.)','FontSize',25);xlabel('Depth(m)','FontSize',15);
ylabel('Arrival time(ms)','FontSize',15);
grid(gca,'minor');
set(gca,'YDir','reverse');

subplot(1,3,3);
scatter(well_z_cs,well_t_cs_corr*1000,15,'filled','MarkerEdgeColor','b','MarkerFaceColor','c');ylim([0 max(td_pairs(:,2))*1000*1.15]);
title('Check shot(Corrected)','FontSize',25);xlabel('Depth(m)','FontSize',15);
ylabel('Arrival time(ms)','FontSize',15);
grid(gca,'minor');
set(gca,'YDir','reverse');

figure(3);
plot(well_path_vel,well_path_zz,'--r','LineWidth',3); % Plot interval vels seen by well path (according to deviation survey)
hold on;
plot(vel_cs_raw,well_z_cs,'-b','LineWidth',3);  % Plot check shot interval vels
hold on;
plot(vel_cs_corr,well_z_cs,'-g','LineWidth',3);  % Plot check shot interval vels
legend;
hold off;
title('Check shot','FontSize',25);xlabel('Velocity(m/s)','FontSize',15);xlim([200 7500]);%ylim([well_top_z well_bot_z]);
ylabel('Depth(m)','FontSize',15);
grid(gca,'minor');
set(gca,'YDir','reverse');

hold off;


%% --------------------------------------------------------
t2=clock;
deltime=etime(t2,t1);
disp(['raytrace time ' num2str(deltime) ' seconds']);


end
function[sx sy]=linear_eq_solve(xa1,ya1,xa2,ya2,xb1,yb1,xb2,yb2)

ma = (ya2-ya1)/(xa2-xa1);% Slope of line a
mb = (yb2-yb1)/(xb2-xb1);% Slope of line b
ca=(ya1*xa2-ya2*xa1)/(xa2-xa1);% Intercept of line a
cb=(yb1*xb2-yb2*xb1)/(xb2-xb1);% Intercept of line b

sx = (ca-cb)/(mb-ma); % Solution x
sy=(ma*cb-mb*ca)/(ma-mb);% Solution y

end

