function [vol_3d_orig,shifts] = StructuralTensor_3Ddip_cgls

%% Parameters
 
plot_on = 1; % make some plots
sigma = 3; % controls size and shape of Gaussian filter
scale_sigma = 2;

%% Load a 3D volume from OpendTect. Edit odmex.par
vol_3d_orig = readcbvs('/apps/gsc/matlab-library/opendtect_link/odmex.par');
vol_3d_orig = single(vol_3d_orig);
vol_3d_scale = vol_3d_orig;
vol_3d_scale = sign(vol_3d_scale).*log(1+abs(vol_3d_scale));
%vol_3d = permute(vol_3d,[1 3 2]);
[nt,nxl,nil] = size(vol_3d_scale);
nsamp = nt*nxl*nil;

disp('Calculate gradients')
% Calculate gradients in x, y and z directions
[Ix, Iz, Iy] = gradient(vol_3d_scale); % Ix crossline, Iz time, Iy inline 
% The first output FX is always the gradient
% along the 2nd dimension of F, going across columns.
% The second output FY is always the gradient along
% the 1st dimension of F, going across rows.  For
% the third output FZ and the outputs that follow,
% the Nth output is the gradient along the Nth dimension of F.

% Calculate gradients on vector
% vol_3d_v = vol_3d(:);
clear vol_3d_scale

Ixy = Ix.*Iy;
Ixz = Ix.*Iz;
Iyz = Iy.*Iz;
Iyy = Iy.*Iy;
Ixx = Ix.*Ix;
Izz = Iz.*Iz;

clear Ix Iy Iz

% Smooth gradient components with a Gaussian
Ixy = imgaussian(Ixy,sigma,scale_sigma*sigma);
Ixz = imgaussian(Ixz,sigma,scale_sigma*sigma);
Iyz = imgaussian(Iyz,sigma,scale_sigma*sigma);
Iyy = imgaussian(Iyy,sigma,scale_sigma*sigma);
Ixx = imgaussian(Ixx,sigma,scale_sigma*sigma);
Izz = imgaussian(Izz,sigma,scale_sigma*sigma);

Ixy = Ixy(:);
Ixz = Ixz(:);
Iyz = Iyz(:);
Iyy = Iyy(:);
Ixx = Ixx(:);
Izz = Izz(:);

disp('Calculate slopes')
for i_loop = 1:length(Ixx);
    [l1,l2,~,nx(i_loop,1),ny(i_loop,1),nz(i_loop,1),~,~,~,~,~,~] = EigenVectors3D(Ixx(i_loop), Ixy(i_loop), Ixz(i_loop), Iyy(i_loop), Iyz(i_loop), Izz(i_loop));
    % [~,~,~,e1{i_loop},~,~] = EigenVectors3D([Ixx(i_loop), Ixy(i_loop) , Ixz(i_loop); Ixy(i_loop) , Iyy(i_loop) , Iyz(i_loop); Ixz(i_loop) , Iyz(i_loop) , Izz(i_loop)]);   
    % [e1{i_loop},e2{i_loop},e3{i_loop},l1(i_loop),l2(i_loop),l3(i_loop)] = eigen_decomposition([Ixx(i_loop), Ixy(i_loop) , Ixz(i_loop); Ixy(i_loop) , Iyy(i_loop) , Iyz(i_loop); Ixz(i_loop) , Iyz(i_loop) , Izz(i_loop)]);   
    planarity(i_loop,1) = (l1-l2)./l1;
end

% disp('Calculate slopes')
% [nx,ny,nz,planarity] = calculate_dip3d(Ix,Iy,Iz,sigma,scale_sigma);   
clear I*

disp('Calculate shifts')
tol = 1e-4;
maxit = 200;
eta = 0.1;

data = make_data(nsamp,nt,nxl,nil,double(nx),double(ny),double(nz),planarity,eta);

[shifts,~,~,~] = cgls(@(x,flag)calcshifts(x,flag,nsamp,nt,nil,nxl,nx,ny,double(nz),planarity,eta)...
    ,data,0,tol,maxit,'true',[]);
    
shifts = reshape(shifts,nt,nxl,nil);

z_cube = repmat(1:nt,nxl,1)';
z_cube = repmat(z_cube,1,nil);
z_cube = reshape(z_cube,nt,nxl,nil);

shifts = z_cube+shifts;

plot_event(150,vol_3d_orig,z_cube,shifts,nt,nxl,nil)

end

function plot_event(event,vol_3d_orig,z_cube,shifts,nt,nxl,nil)
    figure
    event = z_cube.*(shifts > event & shifts < event+1);

    figure
    subplot(1,3,1); imagesc(squeeze(vol_3d_orig(:,:,5)))
    subplot(1,3,2); imagesc(squeeze(shifts(:,:,5)))
    subplot(1,3,3); imagesc(squeeze(vol_3d_orig(:,:,5))); hold on; h = imagesc(1e5.*ones(nt,nxl)); hold off; set(h, 'AlphaData', event(:,:,5))
end

function y = calcshifts(x,flag,nsamp,nt,nil,nxl,nx,ny,nz,planarity,eta) 
    if flag == 1; % A*x
        y = (spdiags([planarity;planarity;repmat(eta,nsamp,1)],0,3*nsamp,3*nsamp)...             %*spdiags([nz;nz],0,2*nsamp,2*nsamp)...
            *[spdiags(nz,0,nsamp,nsamp)*calculate_dx(nt,nxl,nil,'no_transp');spdiags(nz,0,nsamp,nsamp)*calculate_dy(nt,nxl,nil,'no_transp');calculate_dz(nt,nxl,nil,'no_transp')])...
            *x;    
    elseif flag == 2; % A'*x
        y = (spdiags([planarity;planarity;repmat(eta,nsamp,1)],0,3*nsamp,3*nsamp)... %*spdiags([nz;nz],0,2*nsamp,2*nsamp)...
            *[spdiags(nz,0,nsamp,nsamp)*calculate_dx(nt,nxl,nil,'no_transp');spdiags(nz,0,nsamp,nsamp)*calculate_dy(nt,nxl,nil,'no_transp');calculate_dz(nt,nxl,nil,'no_transp')])'...
            *x;
    end
end

function data = make_data(nsamp,nt,nxl,nil,nx,ny,nz,planarity,eta)     
    data = (spdiags([planarity;planarity;repmat(eta,nsamp,1)],0,3*nsamp,3*nsamp)*[nx;ny;zeros(nsamp,1)]);
end

function dx = calculate_dx(nt,nxl,nil,transp_flag) % x crossline
    %data should be a column vector 
    %from 3d volume with dimensions nt,nxl,nil
   if strcmp(transp_flag,'no_transp')
        % dx
        nsamp = nt*nxl*nil; 
        dx1 = [-1*ones(nt,1); zeros((nxl-2)*nt,1); ones(nt,1)];
        dx2 = [ones(nt,1); 1/2*ones((nxl-2)*nt,1); zeros(nt,1)];
        dx2 = [ones(nt,1); dx2(1:end-nt,1)];        
        dx3 = [-1/2*ones((nxl-2)*nt,1); -1*ones(nt*2,1)];          
       
        dx = spdiags([dx3,dx1,dx2],[-nt 0 nt],nt*nxl,nt*nxl);
        
        for i_loop = 1:nil;
            dx_cell{i_loop} = dx;
        end
                
        dx = blkdiag(dx_cell{:});
        %dx = dx*data;
        
        %s_dx= dx*double(vol_3d_v);
        %s_dx = reshape(s_dx,nt,nxl,nil);        
   elseif strcmp(transp_flag,'transp')
        % dx' 
        nsamp = nt*nxl*nil;
        dx1 = [-1*ones(nt,1); zeros((nxl-2)*nt,1); ones(nt,1)];
        dx2 = [ones(nt,1); 1/2*ones((nxl-2)*nt,1); zeros(nt,1)];
        %dx2 = [ones(nt,1); dx2(1:end-nt,1)];        
        dx3 = [-1/2*ones((nxl-2)*nt,1); -1*ones(nt*2,1)];    
        dx = spdiags([dx2,dx1,dx3],[-nt 0 nt],nt*nxl,nt*nxl);
        
        for i_loop = 1:nil;
            dx_cell{i_loop} = dx;
        end
                
        dx = blkdiag(dx_cell{:});
        %dx = dx*data;
   end
end

function dy = calculate_dy(nt,nxl,nil,transp_flag) % y inline
    % data should be a column vector 
    % from 3d volume with dimensions nt,nxl,nil
    if strcmp(transp_flag,'no_transp')
        nsamp = nt*nxl*nil; 
        dy1 = [-1*ones(nt*nxl,1); zeros(nsamp-2*nt*nxl,1); ones(nt*nxl,1)];
        dy2 = [ones(2*nt*nxl,1); 1/2*ones(nt*nxl*nil-2*nt*nxl,1);];
        dy3 = [-1/2*ones(nt*nxl*nil-2*nt*nxl,1); -1*ones(2*nt*nxl,1)];       
      
        dy = spdiags([dy3,dy1,dy2],[-nxl*nt 0 nxl*nt],nt*nxl*nil,nt*nxl*nil);             
        %dy = dy*data;
    elseif strcmp(transp_flag,'transp')
        nsamp = nt*nxl*nil; 
        dy1 = [-1*ones(nt*nxl,1); zeros(nsamp-2*nt*nxl,1); ones(nt*nxl,1)];
        dy2 = [ones(nt*nxl,1); 1/2*ones(nt*nxl*nil-nt*nxl,1);];
        dy3 = [-1/2*ones(nt*nxl*nil-nt*nxl,1); -1*ones(nt*nxl,1)];         
      
        dy = spdiags([dy2,dy1,dy3],[-nxl*nt 0 nxl*nt],nt*nxl*nil,nt*nxl*nil);             
        %dy = dy*data;    
    end
end

function dz = calculate_dz(nt,nxl,nil,transp_flag)
        % data should be a column vector 
        % from 3d volume with dimensions nt,nxl,nil
        if strcmp(transp_flag,'no_transp')
            % dz
            nsamp = nt*nxl*nil;       
            dz1 = repmat([-1 zeros(1,nt-2) 1],1,nsamp/nt)';
            dz2 = repmat([-1/2*ones(1,nt-2) -1 0],1,nsamp/nt)';
            dz3 = repmat([1 1/2*ones(1,nt-2) 0],1,nsamp/nt)';
            dz3 = [1;dz3(1:end-1)];
            dz = spdiags([dz2,dz1,dz3],[-1 0 1],nsamp,nsamp);
            %dz = dz*data;
        elseif strcmp(transp_flag,'transp')
            % dz'
            nsamp = nt*nxl*nil;
            dz1 = repmat([-1 zeros(1,nt-2) 1],1,nsamp/nt)';
            dz2 = repmat([-1/2*ones(1,nt-2) -1 0],1,nsamp/nt)';
            dz3 = repmat([1 1/2*ones(1,nt-2) 0],1,nsamp/nt)';       
            dz2 =  [-1/2;dz2(1:end-1)];
            dz = spdiags([dz3,dz1,dz2],[-1 0 1],nsamp,nsamp);
            %dz = dz*data;
        end      
end



