function []=ray_trace_vel_model(volume_path,il_byte_loc,xl_byte_loc,x_interval,direction)
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
s_a1 = 8 ;             % to control the spread of the ray fan. A bigger number means bigger spread
s_a2= 1000;             % to control the decimation of the ray fan
r_vel = 1540;           % Mimimum rock velocity
% Plot parameters
offset = 9000;         % Maximum Offset in Plots
z_max=4500;             % Maximum Z in Ray Tracing Plot
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
[c,i] = min(v_temp(5:size(v_temp,1),:));
z_wb=z(i);
clear v_temp c;
%---------------------------------------------------------------------------------------
nx=length(x);
nz=length(z);
dx=x(2)-x(1);
dg=z(2)-z(1);

rayvelmod_bg(vel,dx,dg);                        % Install the velocity model
close all;
figure(2);
imagesc(x,z,vel,[1400 v_max]);ylim([0 z_max]);   %Plot the velocity model
hold all;
%plot(x,z_wb,'-k','LineWidth',2);                % Plot water bottom

%% ------------------------------------------------------------
% disp(' ');
% disp(' ');
% plotimage(vel-mean(vel(:)),z,x);
% xlabel('meters');ylabel('meters');
% disp(' ');
% disp(' ');
% disp('Consider this velocity model');
% disp(['Solid black is ' int2str(round(max(vel(:)))) ' m/s']);
% disp(['Solid white is ' int2str(round(min(vel(:)))) ' m/s']);
x0=round(nx/3)*dx;z0=0;

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


%------------------Estimate tmax,dt,tstep----------------------
vlow=min(min(vel));
tmax=s*max(z)/vlow;dt=.004;tstep=0:dt:tmax;
angles=(-pi/2+pi/s_a1):pi/s_a2:(pi/2-pi/s_a1);
%-----------------Ray Tracing------------------------------------------------
indx=near(x,x0);indz=near(z,z0);
v0=vel(indz,indx);
sr=zeros(length(angles),5);% Matrix to keep shot record
%trace the rays
t1=clock;
k_min=20;
for k=1:length(angles)
    r0=[x0 z0 sin(angles(k))/v0 cos(angles(k))/v0];     % Define the initial ray vector
    [t,r]=shootrayvxz_bg(tstep,r0);                     % Shoot the rays
    line(r(:,1),r(:,2),ones(size(t)),'color','w');      % Plot the rays
    sr(k,1)=angles(k);                                  % Store angles a first column
    z_k=r(:,2);                                         % This is when the z coordinates in the ray's travel path 
    [z_min ind1]=min(abs(z_k(k_min:length(z_k))));      % %remove first couple of points near datum and Find when the ray reaches the datum
    
    % If the ray reaches the datum at all
    if (z_min<10)
        sr(k,2)=r((k_min+ind1-1),1);                    % Store the x location at which the ray is received
        sr(k,3)=t(k_min+ind1-1);                        % Store the the total travel time when the ray is received
        [sr(k,4),i]=max(z_k(1:length(z_k)));            % Store the the maximum z the ray sees
        [c,wb_x_index]=min(x-r(i,1));                   % Index of x-axis at the location of maximum depth reached
        wb_z_k=z_wb(wb_x_index);                        % Water bottom depth correspoinding to maximum depth reached
        sr(k,5)=sr(k,4)-wb_z_k;                         % the overburden thickness correspongin to the maximum depth reached
    end
end

clear c wb_z_k wb_x_index z_min ind1;
title('Ray Tracing','FontSize',25); xlabel('x axis (m)','FontSize',15); ylabel('z axis (m)','FontSize',15);
set(gca,'XAxisLocation','top');ylim([0 z_max]);

%------------------Plot results------------------------------------
figure(3);

subplot(1,3,1);
scatter(sr(:,2)-x0,sr(:,3)*1000,20,'filled','MarkerEdgeColor','b','MarkerFaceColor','c');
title('Shot record (z)','FontSize',25); xlabel('Offset (m)','FontSize',15); 
ylabel('Arrival time(ms)','FontSize',15);xlim([-offset offset]);ylim([2000 6500]);
grid(gca,'minor');
set(gca,'YDir','reverse');

subplot(1,3,2);
scatter(sr(:,2)-x0,sr(:,4),20,'filled','MarkerEdgeColor','r','MarkerFaceColor','y');
title('Max Depth Reached','FontSize',25); xlabel('Offset (m)','FontSize',15); 
ylabel('Maximum Depth)','FontSize',15);xlim([-offset offset]);ylim([1000 4500]);
grid(gca,'minor');
set(gca,'YDir','reverse');

subplot(1,3,3);
scatter(sr(:,2)-x0,sr(:,5),20,'filled','MarkerEdgeColor','k','MarkerFaceColor','m');
title('Max Depth Overburden Reached','FontSize',25); xlabel('Offset (m)','FontSize',15); 
ylabel('Maximum Depth)','FontSize',15);xlim([-offset offset]);ylim([-100 4000]);
grid(gca,'minor');
set(gca,'YDir','reverse');
%-------------------------------------------------------------------------

t2=clock;
deltime=etime(t2,t1);
disp(['raytrace time ' num2str(deltime) ' seconds']);


end
