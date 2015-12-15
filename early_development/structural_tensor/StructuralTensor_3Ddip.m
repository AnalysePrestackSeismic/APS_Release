function StructuralTensor_3Ddip

%% Parameters
 
plot_on = 1; % make some plots
sigma = 5; % controls size and shape of Gaussian filter
scale_sigma = 4;

%% Load a 3D volume from OpendTect. Edit odmex.par
%vol_3d_orig = readcbvs('/bgdata/matlab_opendtect/mexfiles/odmex.par');
vol_3d_orig = readcbvs('/apps/gsc/matlab-library/opendtect_link/odmex.par');
vol_3d_orig = single(vol_3d_orig);
vol_3d = vol_3d_orig;
vol_3d = sign(vol_3d).*log(1+abs(vol_3d));
%vol_3d = permute(vol_3d,[1 3 2]);
[nt,nxl,nil] = size(vol_3d);

%imagesc(reshape(vol_3d_v(1:nt*nxl),nt,nxl))

if plot_on == 1
    slice_no = floor(nt/2);
    figure(1)
    subplot(3,2,1); imagesc(squeeze(vol_3d(slice_no,:,:))); title(['Data - Slice no. ',num2str(slice_no)]); 
end

%% Calculate Structural Tensor

% Smooth the 3D seismic volume
% vol_3d = imgaussian(vol_3d,sigma,scale_sigma*sigma);

if plot_on == 1
    slice_no = floor(nt/2);
    figure(1)
    subplot(3,2,2); imagesc(squeeze(vol_3d(slice_no,:,:))); title(['Smooth data - Slice no. ',num2str(slice_no)]);
end

% Calculate gradients in x, y and z directions
[Ix, Iz, Iy] = gradient(vol_3d); % Ix crossline, Iz time, Iy inline 
% The first output FX is always the gradient
% along the 2nd dimension of F, going across columns.
% The second output FY is always the gradient along
% the 1st dimension of F, going across rows.  For
% the third output FZ and the outputs that follow,
% the Nth output is the gradient along the Nth dimension of F.

% Calculate gradients on vector
vol_3d_v = vol_3d(:);

if plot_on == 1
    slice_no = floor(nt/2);
    figure(1)
    subplot(3,2,3); imagesc(squeeze(Ix(slice_no,:,:))); title(['Ix - Slice no. ',num2str(slice_no)]);
    subplot(3,2,4); imagesc(squeeze(Iy(slice_no,:,:))); title(['Iy - Slice no. ',num2str(slice_no)]);
    subplot(3,2,5); imagesc(squeeze(Iz(slice_no,:,:))); title(['Iz - Slice no. ',num2str(slice_no)]);
end

[nx,ny,nz,planarity] = calculate_dip3d(Ix,Iy,Iz,sigma,scale_sigma);   
    % [e1{i_loop},e2{i_loop},e3{i_loop},l1(i_loop),l2(i_loop),l3(i_loop)] = eigen_decomposition([Ixx(i_loop), Ixy(i_loop) , Ixz(i_loop); Ixy(i_loop) , Iyy(i_loop) , Iyz(i_loop); Ixz(i_loop) , Iyz(i_loop) , Izz(i_loop)]);   

%p_vol = reshape(p,nt,nil,nxl);
%q_vol = reshape(q,nt,nil,nxl);

% p = -nx/nt (inline slope)
% q = -ny/nt (crossline slope)
%slopes = [p(:); q(:)];

% Solve using PCG
% PCG requires A'A as input
data = double([nx(:);ny(:)]);
nz = double(nz);
planarity = double(planarity);

tol = 1e-6;
maxit = 1000;
eta = 0.1;
x0 = zeros(nt*nil*nxl,1);

[shifts] = pcg_edit(@(x)calcshifts(x,nz,nt,nil,nxl,planarity,eta),[planarity(:);planarity(:)].*data,tol,maxit,[],[],x0);

shifts = reshape(shifts,nt,nxl,nil);

% data = double([nx(:)./nz(:);ny(:)./nz(:)]);
% nz = double(nz)';
% planarity = double(planarity);
% 
% tol = 1e-6;
% maxit = 1000;
% eta = 0.1;
% x0 = zeros(nt*nil*nxl,1);
% 
% dx = calculate_dx(nt,nxl,nil,'no_transp');
% dy = calculate_dy(nt,nxl,nil,'no_transp');
% %dz = calculate_dz(nt,nxl,nil,'no_transp')
% 
% %dz = eta.*dz;
% 
% A = [dx;dy];
% 
% shifts = lsqr(A,data,tol,maxit,[],[],x0);
% shifts = reshape(shifts,nt,nxl,nil);

z_cube = repmat(1:nt,nxl,1)';
z_cube = repmat(z_cube,1,nil);
z_cube = reshape(z_cube,nt,nxl,nil);

shifts1 = z_cube+shifts;
shifts2 = z_cube-shifts;

figure
subplot(2,1,1); imagesc(squeeze(vol_3d_orig(:,10,:)))
subplot(2,1,2); imagesc(squeeze(shifts1(:,10,:)))

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

function [y] = calcshifts(x,nz,nt,nil,nxl,planarity,eta)
    % Calculate the result of 
    dx = calculate_dx(nt,nxl,nil,'no_transp');
    dy = calculate_dy(nt,nxl,nil,'no_transp');
    
    % Calculate gradient
    y = spdiags([planarity;planarity],0,2*nt*nxl*nil,2*nt*nxl*nil)*([nz;nz].*([dx;dy]*x));
end

function [g] = tensfiltA(data,nX,nZ,u2)
% apply A'A filter from D. Hale, 2007, CWP Report 567, "Local dip filtering
% with directional Laplacians"
    for ii = 2:nZ
        for jj=2:nX
            u2i = u2(ii,jj);
            u1i = sqrt(1-u2i*u2i);
            a11 = 1-u1i*u1i;
            a12 = -u1i*u2i;
            a22 = 1-u2i*u2i;
            fa = data(ii,jj)-data(ii-1,jj-1);
            fb = data(ii,jj-1)-data(ii-1,jj);
            f1 = 0.5*(fa-fb);
            f2 = 0.5*(fa+fb);
            g1 = a11*f1+a12*f2;
            g2 = a12*f1+a22*f2;
            ga = 0.5*(g1+g2);
            gb = 0.5*(g1-g2);
            g(ii,jj) = ga;
            g(ii-1,jj-1) = g(ii-1,jj-1)-ga;
            g(ii,jj-1) = g(ii,jj-1)-gb;
            g(ii-1,jj) = g(ii-1,jj)+gb;
        end
    end
end

function [y] = smooth(x,alpha,nX,nZ,u2)
% smoothing by inversion regularization using Claerbout's wavekill filter
% (filter A from D. Hale, 2007, CWP Report 567, "Local dip filtering
% with directional Laplacians"). Alpha controls smoothness.
    x = reshape(x,nZ,nX);
    y = reshape(x+alpha.*tensfiltA(x,nX,nZ,u2),[],1);
end



