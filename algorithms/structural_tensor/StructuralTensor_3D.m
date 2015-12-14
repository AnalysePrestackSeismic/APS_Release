function StructuralTensor_3D

%% Parameters
 
plot_on = 1; % make some plots
sigma = 2; % controls size and shape of Gaussian filter
scale_sigma = 4;
matlab_load = 1;
% if matlab_load == 1;
%     %%
%     job_meta = load(job_meta_path);                                 % Load job meta file
% 
%     pkey_inc_mode = mode(job_meta.pkey_inc);                        % Primary Key (inline) increment (mode )
%     skey_inc_mode = mode(job_meta.skey_inc);                        % Secondary Key (inline) Increment (mode)
% 
%     pkeyn = 1+((job_meta.pkey_max(str2double(vol_index))-job_meta.pkey_min(str2double(vol_index)))...
%         /job_meta.pkey_inc(str2double(vol_index)));                 % Calculate Number of inlines (primary key)
%     skeyn = 1+((job_meta.skey_max(str2double(vol_index))-job_meta.skey_min(str2double(vol_index)))...
%         /job_meta.skey_inc(str2double(vol_index)));                 % Calculate Number of inlines (secondary key)
% 
%     if str2double(end_slab) > job_meta.n_samples{str2double(vol_index)}
%         end_slab = num2str(job_meta.n_samples{str2double(vol_index)});                    % update endslab if the data provided is of shorter length
%     end
% 
%     n_slices = str2double(end_slab)-str2double(start_slab)+1;   % Number of Slices
%     vol_3d = zeros(pkeyn*skeyn,n_slices,'single');                                      % Initalize matrix for all inlines, xlines and slices
% 
%     loopfin = size(job_meta.liveblocks,1);                          % Number of live blocks
%     lpi = 1;
% 
%     while lpi <= loopfin
%         i_block = job_meta.liveblocks(lpi);                            % Block Number for Current Live Block
%         [~, traces, ilxl_read, ~] = ...
%             node_segy_read(job_meta_path,vol_index,num2str(i_block));
%         traces = [zeros(1,size(traces,2)); traces(2:end,:)];
% 
%         n_iline = (ilxl_read(:,1)-job_meta.pkey_min(str2double(vol_index)))/pkey_inc_mode+1;
%         n_xline = (ilxl_read(:,2)-job_meta.skey_min(str2double(vol_index)))/skey_inc_mode+1;
%         lin_ind = ((n_iline-1).*skeyn)+n_xline;    
% 
%         %slices = zeros(pkeyn*skeyn,n_slices,'single')
%         vol_3d(double(lin_ind),1:n_slices) = traces(str2double(start_slab):str2double(end_slab),:)';        
% 
%         fprintf('-- Block %d of %d --\n',lpi,loopfin);
%         lpi = lpi + 1;
% 
%         % Save water bottom pick    
%     end
%     
%     vol_3d = reshape(vol_3d,skeyn,pkeyn,n_slices);
%     vol_3d = single(vol_3d);    
% else
    %% Load a 3D volume from OpendTect. Edit odmex.par
    vol_3d = readcbvs('/apps/gsc/matlab-library/opendtect_link/odmex.par');
    vol_3d = single(vol_3d);
%end

vol_3d = sign(vol_3d).*log(1+abs(vol_3d));

[n_samp,n_xline,n_iline] = size(vol_3d);

if plot_on == 1
    slice_no = floor(n_samp/2);
    figure(1)
    subplot(3,2,1); imagesc(squeeze(vol_3d(slice_no,:,:))); title(['Data - Slice no. ',num2str(slice_no)]);
end

%% Calculate Structural Tensor

% Smooth the 3D seismic volume
% vol_3d = imgaussian(vol_3d,sigma,scale_sigma*sigma);
if plot_on == 1
    slice_no = floor(n_samp/2);
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

% Ix = imgaussian(Ix,sigma,scale_sigma*sigma);
% Iy = imgaussian(Iy,sigma,scale_sigma*sigma);
% Iz = imgaussian(Iz,sigma,scale_sigma*sigma);

%clear vol_3d

if plot_on == 1
    slice_no = floor(n_samp/2);
    figure(1)
    subplot(3,2,3); imagesc(squeeze(Ix(slice_no,:,:))); title(['Ix - Slice no. ',num2str(slice_no)]);
    subplot(3,2,4); imagesc(squeeze(Iy(slice_no,:,:))); title(['Iy - Slice no. ',num2str(slice_no)]);
    subplot(3,2,5); imagesc(squeeze(Iz(slice_no,:,:))); title(['Iz - Slice no. ',num2str(slice_no)]);
end

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

% if ~exist('parpool') 
%     parpool
% end

% Calculate eigenvectors and eigenvalues

%matlabpool('local')


% parfor i_loop = 1:length(Ixx);
%     [e1{i_loop},e2{i_loop},e3{i_loop},l1(i_loop),l2(i_loop),l3(i_loop)] = eigen_decomposition([Ixx(i_loop), Ixy(i_loop) , Ixz(i_loop); Ixy(i_loop) , Iyy(i_loop) , Iyz(i_loop); Ixz(i_loop) , Iyz(i_loop) , Izz(i_loop)]);   
%     % [e1{i_loop},e2{i_loop},e3{i_loop},l1(i_loop),l2(i_loop),l3(i_loop)] = eigen_decomposition([Ixx(i_loop), Ixy(i_loop) , Ixz(i_loop); Ixy(i_loop) , Iyy(i_loop) , Iyz(i_loop); Ixz(i_loop) , Iyz(i_loop) , Izz(i_loop)]);   
% end
%parpool

% parfor i_loop = 1:length(Ixx);
%     [x_dip(i_loop),y_dip(i_loop)] = calculate_dip([Ixx(i_loop), Ixy(i_loop), Ixz(i_loop); Ixy(i_loop), Iyy(i_loop), Iyz(i_loop); Ixz(i_loop), Iyz(i_loop), Izz(i_loop)]);   
%     % [e1{i_loop},e2{i_loop},e3{i_loop},l1(i_loop),l2(i_loop),l3(i_loop)] = eigen_decomposition([Ixx(i_loop), Ixy(i_loop) , Ixz(i_loop); Ixy(i_loop) , Iyy(i_loop) , Iyz(i_loop); Ixz(i_loop) , Iyz(i_loop) , Izz(i_loop)]);   
% end

parfor i_loop = 1:length(Ixx);
    [~,~,~,e1x(i_loop),e1y(i_loop),e1z(i_loop),~,~,~,~,~,~] = EigenVectors3D(Ixx(i_loop), Ixy(i_loop), Ixz(i_loop), Iyy(i_loop), Iyz(i_loop), Izz(i_loop));
    % [~,~,~,e1{i_loop},~,~] = EigenVectors3D([Ixx(i_loop), Ixy(i_loop) , Ixz(i_loop); Ixy(i_loop) , Iyy(i_loop) , Iyz(i_loop); Ixz(i_loop) , Iyz(i_loop) , Izz(i_loop)]);   
    % [e1{i_loop},e2{i_loop},e3{i_loop},l1(i_loop),l2(i_loop),l3(i_loop)] = eigen_decomposition([Ixx(i_loop), Ixy(i_loop) , Ixz(i_loop); Ixy(i_loop) , Iyy(i_loop) , Iyz(i_loop); Ixz(i_loop) , Iyz(i_loop) , Izz(i_loop)]);   
end

clear I*

e1x = reshape(e1x,n_samp,n_xline,n_iline);
e1y = reshape(e1y,n_samp,n_xline,n_iline);
e1z = reshape(e1z,n_samp,n_xline,n_iline);

x_dip = atand(e1x./e1z);
y_dip = atand(e1y./e1z);

fid = fopen('/data/TZA/dtect/2014_TZA_block_34_mafia/Misc/x_dip_2.bin','w');
fwrite(fid,x_dip,'float32');
fclose(fid);
fid = fopen('/data/TZA/dtect/2014_TZA_block_34_mafia/Misc/y_dip_2.bin','w');
fwrite(fid,y_dip,'float32');
fclose(fid);
figure

% subplot(3,2,1); imagesc(squeeze(atand(e1z(:,:,10)./squeeze(e1x(:,:,10))))-squeeze(atand(e2z(:,:,10)./squeeze(e2x(:,:,10)))))
% subplot(3,2,2); imagesc(squeeze(sign(e1x(:,:,10))))
% subplot(3,2,3); imagesc(squeeze(atand(e1z(:,:,10)./squeeze(e1x(:,:,10)))))
% subplot(3,2,4); imagesc(squeeze(atand(e2z(:,:,10)./squeeze(e2x(:,:,10)))),[-45 45])
% subplot(3,2,5); imagesc(squeeze(atand(e1z(:,:,10)./squeeze(e1x(:,:,10))))-(squeeze(sign(e1x(:,:,10))).*90),[-45 45])
% subplot(3,2,6); imagesc(squeeze(atand(e3z(:,:,10)./squeeze(e3x(:,:,10)))),[-45 45])

% l1 = reshape(l1,n_samp,n_xline,n_iline);
% l2 = reshape(l2,n_samp,n_xline,n_iline);
% l3 = reshape(l3,n_samp,n_xline,n_iline);

%% Attributes
% C_plane = (l1-l2)./(l1+l2);
% C_line = (l2-l3)./(l2+l3);
% C_fault = (2*l2.*(l2-l3))./((l1+l2).*(l2+l3));

% e1 = cell2mat(e1);
% e1x = reshape(e1(1,:),n_samp,n_xline,n_iline);
% e1y = reshape(e1(2,:),n_samp,n_xline,n_iline);
% e1z = reshape(e1(3,:),n_samp,n_xline,n_iline);
% 
% e1 = cell2mat(e2);
% clear e2
% e2x = reshape(e1(1,:),n_samp,n_xline,n_iline);
% e2y = reshape(e1(2,:),n_samp,n_xline,n_iline);
% e2z = reshape(e1(3,:),n_samp,n_xline,n_iline);
% 
% e1 = cell2mat(e3);
% clear e3
% e3x = reshape(e1(1,:),n_samp,n_xline,n_iline);
% e3y = reshape(e1(2,:),n_samp,n_xline,n_iline);
% e3z = reshape(e1(3,:),n_samp,n_xline,n_iline);
% clear e1

%% Diffusion tensor filtering

% e1{i_loop} - largest gradient 
% e2{i_loop} -
% e3{i_loop} - 

% [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz] = Structure2Diffusion(e1,e2,e3,l1,l2,l3,n_samp,n_xline,n_iline);
% 
% T = 2;
% dt = 0.15;
% dt_max = dt;
% t = 0;
% while (t < (T-0.001))
%     dt = min(dt_max,T-t); 
%     t = t + dt;
%  
%     vol_3d_sm = diffusion_scheme_3D(vol_3d_orig,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz,dt);
% 
% end
% 
% subplot(3,1,1); imagesc(squeeze(vol_3d_orig(:,:,10)));
% subplot(3,1,2); imagesc(squeeze(vol_3d_sm(:,:,10)));
% subplot(3,1,3); imagesc(squeeze(vol_3d_orig(:,:,10))-squeeze(vol_3d_sm(:,:,10)));
% 
% y1 = pcg(@(x)smooth(x,3,nX,nZ,sign(tens(:,:,2)).*tens(:,:,3)),data(:),[],500);
% 
% @(x)smooth(vol_3d_orig,3,nX,nZ,sign(e1x).*e1z)

% C_fault = (2*l2.*(l2-l3))./((l1+l2).*(l2+l3));
% figure(2)
% subplot(2,1,1); imagesc(squeeze(vol_3d(:,:,10)));
% subplot(2,1,2); imagesc(squeeze(C_fault(:,:,10)));

% subplot(3,2,1); imagesc(squeeze(atand(e1z(:,:,10)./squeeze(e1x(:,:,10))))-squeeze(atand(e2z(:,:,10)./squeeze(e2x(:,:,10)))))
% subplot(3,2,2); imagesc(squeeze(sign(e1x(:,:,10))))
% subplot(3,2,3); imagesc(squeeze(atand(e1z(:,:,10)./squeeze(e1x(:,:,10)))))
% subplot(3,2,4); imagesc(squeeze(atand(e2z(:,:,10)./squeeze(e2x(:,:,10)))),[-45 45])
% subplot(3,2,5); imagesc(squeeze(atand(e1z(:,:,10)./squeeze(e1x(:,:,10))))-(squeeze(sign(e1x(:,:,10))).*90),[-45 45])
% subplot(3,2,6); imagesc(squeeze(atand(e3z(:,:,10)./squeeze(e3x(:,:,10)))),[-45 45])

% quiverplot(squeeze(vol_3d(:,:,10)),squeeze(e2x(:,:,10)).*squeeze(C_plane(:,:,10)),squeeze(e2z(:,:,10)).*squeeze(C_plane(:,:,10)),5,2)

% [g] = tensfiltA(squeeze(vol_3d(:,:,10)),nX,nZ,u2)
end

% function [y] = smooth(x,alpha,nX,nZ,u2)
% % smoothing by inversion regularization using Claerbout's wavekill filter
% % (filter A from D. Hale, 2007, CWP Report 567, "Local dip filtering
% % with directional Laplacians"). Alpha controls smoothness.
%     x = reshape(x,nZ,nX);
%     y = reshape(x+alpha.*tensfiltA(x,nX,nZ,u2),[],1);
% end
% 
% function [y] = sim(data,nX,nZ,tens) 
% % structure-orientated similarity from D. Hale, 2009, CWP Report 635,
% % "Structure-oriented smoothing and semblance").
%     y1 = pcg(@(x)smooth(x,3,nX,nZ,sign(tens(:,:,2)).*tens(:,:,3)),data(:),[],500);
%     y1 = reshape(pcg(@(x)smooth(x,45,nX,nZ,sign(tens(:,:,5)).*tens(:,:,6)),y1.^2,[],500),nZ,nX); 
%     y2 = pcg(@(x)smooth(x,3,nX,nZ,sign(tens(:,:,2)).*tens(:,:,3)),data(:).^2,[],500);
%     y2 = reshape(pcg(@(x)smooth(x,45,nX,nZ,sign(tens(:,:,5)).*tens(:,:,6)),y2,[],500),nZ,nX);
%     y = y1./y2;
%     y(y<0) = 0;
%     y(y>1) = 1;
% end