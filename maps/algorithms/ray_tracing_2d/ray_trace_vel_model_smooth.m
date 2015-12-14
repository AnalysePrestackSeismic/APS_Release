function []=ray_trace_vel_model(volume_path,il_byte_loc,xl_byte_loc,x_interval,direction,output_dir)
%   Inputs:
%         volume_path = Path of seismic Volume
%         il_byte_loc = Inline byte location
%         xl_byte_loc = xlbyte location
%         shot_loc = shote location in [il xl]
%         x_interval = inline or xline interval in m (depending on direction) 
%         direction = 1 if its an inline  0 if its a xline
%         output_dir = output directory
%         smo = smoothening parameter
%%
close all;
[seis_meta il_xl_pos seismic]=segy_to_mat(il_byte_loc,xl_byte_loc,volume_path);

z = 0:(size(seismic,1)-1);
vel =seismic;

if direction==1
    x= double((il_xl_pos (:,2)-il_xl_pos (1,2)))*x_interval;
    %shot_loc = floor( shot_loc(2)-il_xl_pos (1,2))*x_interval;
elseif direction==0
    x= double(il_xl_pos (il_xl_pos (:,1)-il_xl_pos (1,1)))*x_interval;
    %shot_loc = floor( shot_loc(2)-il_xl_pos (1,1))*x_interval;
else
    fprintf ('\nrespecify right dir: 1 if its an inline  0 if its a xline\n')
end
x=x';
z=z'*seis_meta.s_rate/1000;

mat_file_path=strcat(output_dir,'vel_model_mat.mat');



%% ------------------------------------------------------------
%Demo shootrayvxz on the Marmousi model.

%this code finds the file marmousi_mod.mat which should be in the same
%folder as this script
% s=which('raymarmousi_demo');
% ind = findstr(s,'raymarmousi_demo');
% sm=[s(1:ind-1) 'marmousi_mod'];
% disp(['Marmousi model loaded from ' sm])
% load(sm)
nx=length(x);
nz=length(z);
dg=z(2)-z(1);
dx=x(2)-x(1);
pt1=[min(x),min(z)];
pt12=[mean(x),min(z)];
pt2=[max(x),min(z)];
pt23=[max(x),mean(z)];
pt3=[max(x),max(z)];
pt34=[mean(x),max(z)];
pt4=[min(x),max(z)];
pt41=[min(x),mean(z)];

disp(' ')
disp(' ')
plotimage(vel-mean(vel(:)),z,x)
xlabel('meters');ylabel('meters')
disp(' ')
disp(' ')
disp(['Consider this velocity model'])
disp(['Solid black is ' int2str(round(max(vel(:)))) ' m/s'])
disp(['Solid white is ' int2str(round(min(vel(:)))) ' m/s'])

msg='Enter smoother length(meters) (0<=smoother) or -1 or <cr> to end-> :';

r=input(msg);
if(isempty(r)) 
    r=-1; 
end

x0=round(nx/3)*dg;z0=0;

while(r>=0)
    
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

    smooth=r;

    % Define smoother
    nsmooth=2*round(.5*smooth/dg)+1;%odd number
    xb=(0:nx+nsmooth-2)*dx;zb=(0:nz+nsmooth-2)*dg;
    x=(0:nx-1)*dx;z=(0:nz-1)*dg;
    ixcenter=1+(nsmooth-1)/2:nx+(nsmooth-1)/2;
    izcenter=1+(nsmooth-1)/2:nz+(nsmooth-1)/2;

    % run a smoother over it
    t1=clock;
    v=conv2(vel,ones(nsmooth,nsmooth)/(nsmooth*nsmooth),'same');
    t2=clock;
    deltime=etime(t2,t1);
    disp(['smoothing time ' num2str(deltime) ' seconds']);
    plotimage(v-mean(v(:)),z,x)
    xlabel('meters');ylabel('meters')
    title(['Raytracing after ' int2str(smooth) ' m smoother'])
    %install the velocity model
    rayvelmod(v,dg);

    %estimate tmax,dt,tstep
    vlow=min(min(v));
    tmax=max(z)/vlow;dt=.004;tstep=0:dt:tmax;%###############################################?????

    %specify a fan of rays
    %shoot rays towards the side of the model furthest from (x0,z0)
%     d1=norm(pt12-[x0,z0]);
%     d2=norm(pt23-[x0,z0]);
%     d3=norm(pt34-[x0,z0]);
%     d4=norm(pt41-[x0,z0]);
%     if(d1==max([d1 d2 d3 d4]));
%         %shoot to top
%         th1=atan2(pt41(1)-x0,pt41(2)-z0);
%         th2=atan2(pt23(1)-x0,pt23(2)-z0);
%         angles=[th1:(th2-th1)/50:th2];
%     elseif(d2==max([d1 d2 d3 d4]));
%         %shoot to RHS
%         th1=atan2(pt12(1)-x0,pt12(2)-z0);
%         th2=atan2(pt34(1)-x0,pt34(2)-z0);
%         angles=[th1:(th2-th1)/50:th2];
%     elseif(d3==max([d1 d2 d3 d4]));
%         %shoot to bottom
%         th1=atan2(pt41(1)-x0,pt41(2)-z0);
%         th2=atan2(pt23(1)-x0,pt23(2)-z0);
%         angles=[th1:(th2-th1)/50:th2];
%     else
%         %shoot to LHS
%         th1=atan2(pt12(1)-x0,pt12(2)-z0);
%         th2=atan2(pt34(1)-x0,pt34(2)-z0);
%         angles=[th1:(th2-th1)/50:th2];
%     end
    d1=abs(z0-min(z));
    d2=abs(z0-max(z));
    if(d1>=d2);
        %shoot to top
        th1=atan2(pt41(1)-x0,pt41(2)-z0);
        th2=atan2(pt23(1)-x0,pt23(2)-z0);
        angles=[th1:(th2-th1)/100:th2];
    elseif(d2>d1);
        %shoot to bottom
        th1=atan2(pt41(1)-x0,pt41(2)-z0);
        th2=atan2(pt23(1)-x0,pt23(2)-z0);
        angles=[th1:(th2-th1)/100:th2];
    end
    indx=near(x,x0);indz=near(z,z0);
    v0=v(indz,indx);

    %trace the rays
    t1=clock;
    for k=1:length(angles)
        r0=[x0 z0 sin(angles(k))/v0 cos(angles(k))/v0];
        [t,r]=shootrayvxz(tstep,r0);
        line(r(:,1),r(:,2),ones(size(t)),'color','r');
    end
    t2=clock;
    deltime=etime(t2,t1);
    disp(['raytrace time ' num2str(deltime) ' seconds']);
    r=input(msg);
    if(isempty(r)) r=-1; end

end

disp(' ')
disp('You should look at the source file for this demo')
disp('to see how its done. Also type "help raytrace" for more info.') 


%    

end
