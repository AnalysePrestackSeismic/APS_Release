function [vol_3d_orig,shifts] = StructuralTensor_3Ddip_cgls_adpat_smooth

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

%%
disp('Calculate gradients causal')
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
% clear vol_3d_scale

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

% Smooth gradient components with a Gaussian
alpha = 0.2;
Ixy = sm_exponential(Ixy(:),alpha,'causal');
Ixz = sm_exponential(Ixz(:),alpha,'causal');
Iyz = sm_exponential(Iyz(:),alpha,'causal');
Iyy = sm_exponential(Iyy(:),alpha,'causal');
Ixx = sm_exponential(Ixx(:),alpha,'causal');
Izz = sm_exponential(Izz(:),alpha,'causal');

disp('Calculate slopes causal')
for i_loop = 1:length(Ixx);
    [~,~,~,nx_c(i_loop,1),ny_c(i_loop,1),nz_c(i_loop,1),~,~,~,~,~,~] = EigenVectors3D(Ixx(i_loop), Ixy(i_loop), Ixz(i_loop), Iyy(i_loop), Iyz(i_loop), Izz(i_loop));
    % [~,~,~,e1{i_loop},~,~] = EigenVectors3D([Ixx(i_loop), Ixy(i_loop) , Ixz(i_loop); Ixy(i_loop) , Iyy(i_loop) , Iyz(i_loop); Ixz(i_loop) , Iyz(i_loop) , Izz(i_loop)]);   
    % [e1{i_loop},e2{i_loop},e3{i_loop},l1(i_loop),l2(i_loop),l3(i_loop)] = eigen_decomposition([Ixx(i_loop), Ixy(i_loop) , Ixz(i_loop); Ixy(i_loop) , Iyy(i_loop) , Iyz(i_loop); Ixz(i_loop) , Iyz(i_loop) , Izz(i_loop)]);   
end
clear I*

%%
disp('Calculate gradients anticausal')
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
%clear vol_3d_scale

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

Ixy = sm_exponential(Ixy(:),alpha,'anticausal');
Ixz = sm_exponential(Ixz(:),alpha,'anticausal');
Iyz = sm_exponential(Iyz(:),alpha,'anticausal');
Iyy = sm_exponential(Iyy(:),alpha,'anticausal');
Ixx = sm_exponential(Ixx(:),alpha,'anticausal');
Izz = sm_exponential(Izz(:),alpha,'anticausal');

disp('Calculate slopes anticausal')
for i_loop = 1:length(Ixx);
    [~,~,~,nx_a(i_loop,1),ny_a(i_loop,1),nz_a(i_loop,1),~,~,~,~,~,~] = EigenVectors3D(Ixx(i_loop), Ixy(i_loop), Ixz(i_loop), Iyy(i_loop), Iyz(i_loop), Izz(i_loop));
    % [~,~,~,e1{i_loop},~,~] = EigenVectors3D([Ixx(i_loop), Ixy(i_loop) , Ixz(i_loop); Ixy(i_loop) , Iyy(i_loop) , Iyz(i_loop); Ixz(i_loop) , Iyz(i_loop) , Izz(i_loop)]);   
    % [e1{i_loop},e2{i_loop},e3{i_loop},l1(i_loop),l2(i_loop),l3(i_loop)] = eigen_decomposition([Ixx(i_loop), Ixy(i_loop) , Ixz(i_loop); Ixy(i_loop) , Iyy(i_loop) , Iyz(i_loop); Ixz(i_loop) , Iyz(i_loop) , Izz(i_loop)]);   
end

% disp('Calculate slopes')
% [nx,ny,nz,planarity] = calculate_dip3d(Ix,Iy,Iz,sigma,scale_sigma);   
clear I*

%%
uncon_x = reshape((nx_c./nz_c)-(nx_a./nz_a),nt,nxl,nil);
uncon_y = reshape((ny_c./nz_c)-(ny_a./nz_a),nt,nxl,nil);

uncon = uncon_y.*uncon_x;

figure
subplot(1,3,1); imagesc(squeeze(vol_3d_orig(:,:,20)));
subplot(1,3,2); imagesc(squeeze(uncon(:,:,20)),[-0.1 0.1]);
subplot(1,3,3); imagesc(squeeze(vol_3d_orig(:,:,20))); hold on; h = imagesc(1e3.*ones(nt,nxl)); hold off; set(h, 'AlphaData', abs(uncon(:,:,20)) > 0.01)

figure
subplot(1,3,1); imagesc(squeeze(vol_3d_orig(:,50,:)));
subplot(1,3,2); imagesc(squeeze(uncon(:,50,:)));
subplot(1,3,3); imagesc(squeeze(vol_3d_orig(:,50,:))); hold on; h = imagesc(1e3.*ones(nt,nil)); hold off; set(h, 'AlphaData', uncon(:,50,:))

figure
subplot(1,3,1); imagesc(squeeze(vol_3d_orig(150:end,300,:)));
subplot(1,3,2); imagesc(squeeze(uncon(150:end,300,:)),[-0.1 0.1]);

end

function output = sm_exponential(input,a,type)
    % alpha
    b = 1-a;
    
    if strcmp(type,'causal');
        output = [input(1);zeros(numel(input)-1,1)];
        for ii = 2:1:numel(input)
            output(ii) = a*input(ii-1)+b*output(ii-1);
        end
    else strcmp(type,'anticausal');
        output = [zeros(numel(input)-1,1);input(end)]; 
        for ii = numel(input):-1:2
            output(ii) = a*input(ii-1)+b*output(ii-1);
        end
    end
end



