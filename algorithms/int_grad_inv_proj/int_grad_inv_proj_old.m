function [] = int_grad_inv_proj(job_meta_path,i_block)
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
% -------------------------------------------------------------------------
% INT_GRAD_INV_PROJ: function to run Intercept Gradient Inversion using
% dynamic wavelet set.
%   Inputs:
%       seismic_mat_path = path of metadata .mat file.
%       i_block = current block to be processed.
%       n_blocks = total number of blocks to submit for processing.
%   Outputs:
%       digi_intercept = n_blocks of SEGY files.
%       digi_gradient = n_blocks of SEGY files.
%       digi_minimum_energy_eer_projection = n_blocks of SEGY files.
% Authors: Jonathan Edgar and James Selvage
% -------------------------------------------------------------------------
%% Parameters - edit here
output_std = 1;
plot_on = 0;
background = 1;
chi_model_type = 'empirical';
% Load job meta information 
job_meta = load(job_meta_path);
% Tikhonov regularisation weight
wsmooth = 1000;
eer_weight = 1; % Weight for EER constraint
iter = 100;
tol = 1e-2;

%%
% Make ouput directories and create meta information

% Read traces for 2/3 max angle stack
vol_index_wb = ceil(job_meta.nvols*0.6667);
[~, traces{vol_index_wb}, ilxl_read{vol_index_wb}] = ...
    node_segy_read(job_meta_path,num2str(vol_index_wb),i_block);
if size(traces{vol_index_wb},1) == 1 && size(traces{vol_index_wb},2) == 1 
    return
end

% pkey_inc_mode = mode(job_meta.pkey_inc);
% skey_inc_mode = mode(job_meta.skey_inc);
% 
% n_iline = (ilxl_read{vol_index_wb}(:,1)-min(ilxl_read{vol_index_wb}(:,1)))/pkey_inc_mode+1;
% n_xline = (ilxl_read{vol_index_wb}(:,2)-min(ilxl_read{vol_index_wb}(:,2)))/skey_inc_mode+1;
% skeyn = (max(ilxl_read{vol_index_wb}(:,2))-min(ilxl_read{vol_index_wb}(:,2)))/skey_inc_mode+1;
% pkeyn = (max(ilxl_read{vol_index_wb}(:,1))-min(ilxl_read{vol_index_wb}(:,1)))/pkey_inc_mode+1;
% lin_ind = ((n_iline-1).*skeyn)+n_xline;

% Pick water bottom
if isfield(job_meta, 'wb_path')
    %wb_idx = dlmread(job_meta.wb_path,'delimiter','\t');
    wb_idx = dlmread(job_meta.wb_path);
    % col 1 inline
    % col 2 xline
    % col 3 twt
    min_il = min(ilxl_read{vol_index_wb}(:,1));
    max_il = max(ilxl_read{vol_index_wb}(:,1));
    min_xl = min(ilxl_read{vol_index_wb}(:,2));
    max_xl = max(ilxl_read{vol_index_wb}(:,2));
   wb_idx = wb_idx(wb_idx(:,1) >= min_il & wb_idx(:,1) <= max_il & wb_idx(:,2) >= min_xl & wb_idx(:,2) <= max_xl,:); 
%    wb_idx = wb_idx(wb_idx(:,1) >= min_il & wb_idx(:,1) <= (max_il+1) & wb_idx(:,2) >= min_xl & wb_idx(:,2) <= (max_xl+1),:); 

%   wb_idx(:,3) = wb_idx(:,3)./(job_meta.s_rate/1000);
    wb_idx = (wb_idx(:,3)./(job_meta.s_rate/1000))';
    padding = 10;
    wb_idx = wb_idx-padding;    
    win_sub = bsxfun(@plus,wb_idx,(0:job_meta.n_samples{vol_index_wb}-max(wb_idx))');
    win_ind = bsxfun(@plus,win_sub,(0:job_meta.n_samples{vol_index_wb}:...
    job_meta.n_samples{vol_index_wb}*(size(traces{vol_index_wb},2)-1)));
else
    [wb_idx] = water_bottom_picker(traces{vol_index_wb},10);
    % wb_idx_slice = zeros(skeyn,pkeyn);
    % wb_idx_slice(lin_ind) = wb_idx;
    % wb_idx_slice = medfilt3(wb_idx_slice,[3 3]);
    % wb_idx = wb_idx_slice(:);
    wb_idx(wb_idx < 0) = 1;
    win_sub = bsxfun(@plus,wb_idx,(0:job_meta.n_samples{vol_index_wb}-max(wb_idx))');
    win_ind = bsxfun(@plus,win_sub,(0:job_meta.n_samples{vol_index_wb}:...
    job_meta.n_samples{vol_index_wb}*(size(traces{vol_index_wb},2)-1)));
end

% Load block for remaining angle stacks
for i_vol = 1:1:job_meta.nvols 
    if i_vol ~= str2double(vol_index_wb) % don't repeat load for previously read stack
        % Read traces
        [~, traces{i_vol}, ~, ~] = ...
        node_segy_read(job_meta_path,num2str(i_vol),i_block);
    end   
    % Flatten traces to water bottom
    traces{i_vol} = traces{i_vol}(win_ind);    
    %traces{i_vol} = traces{i_vol}(,:);
    input_angles(i_vol) = (job_meta.angle{i_vol}(2)+job_meta.angle{i_vol}(1))/2;
    if plot_on == 1
        figure(1)
        subplot(1,job_meta.nvols,i_vol); imagesc(traces{i_vol});
    end
end

% Test unflatten
% Unflatten data
% [ns,ntraces] = size(traces{1});
% traces_unflat = [traces{1};zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
% for kk = 1:length(wb_idx)
%     traces_unflat(:,kk) = circshift(traces_unflat(:,kk),wb_idx(kk));
% end

%% Load wavelet set
wavelets = load(strcat(job_meta.wav_directory,'all_wavelets_time.mat'));
ns_wavelet = size(wavelets.all_wavelets_time{1},1)-1;
hns_wavelet = floor(ns_wavelet/2);
[ns,ntraces] = size(traces{1});

i_block = str2double(i_block);

for i_vol = 1:1:job_meta.nvols
    wavelet_z_grid = wavelets.all_wavelets_freq{i_vol}(1,:);
    wavelet = wavelets.all_wavelets_freq{i_vol}(2:end,:);
    if plot_on == 1
        figure(2)
        subplot(1,job_meta.nvols,i_vol); imagesc(wavelet);
        figure(3)
        subplot(1,job_meta.nvols,i_vol+1); plot(wavelet(:,10));
    end
    start_interp = min(wavelet_z_grid);%-hns_wavelet;
    end_interp = max(wavelet_z_grid);%+hns_wavelet;
    wavelet_interp{i_vol} = interpolate_wavelets(wavelet_z_grid,wavelet,start_interp,end_interp);
    wavelet_interp{i_vol} = circshift(ifft(wavelet_interp{i_vol}','symmetric'),floor(job_meta.ns_win/2));
    wavelet_interp{i_vol} = wavelet_interp{i_vol};
    
    if start_interp > 1;
        % Pad with zeros or first wavelet
        pad_size = start_interp-1;
        wavelet_interp{i_vol} = [repmat(wavelet_interp{i_vol}(:,1),1,pad_size),...
            wavelet_interp{i_vol}];
    end
    if end_interp < ns;
        % Pad with zeros or last wavelet
        pad_size = ns-end_interp;
        wavelet_interp{i_vol} = [wavelet_interp{i_vol},...
            repmat(wavelet_interp{i_vol}(:,end),1,pad_size)];
    end
    if floor(ns_wavelet/2) == ns_wavelet/2
        wavelet_interp{i_vol}(end+1,:) = wavelet_interp{i_vol}(end,:);
    end
    if plot_on == 1
        figure(4)
        subplot(1,job_meta.nvols,i_vol); imagesc(wavelet_interp{i_vol});
    end
end
interp_wavelet_z_grid = 1:ns;
ns_wavelet = size(wavelet_interp{1},1);
wavelet_norm = wavelet_rms_norm(job_meta,wavelet_interp,interp_wavelet_z_grid,vol_index_wb);
if plot_on == 1
    figure(5)
    imagesc(cell2mat(wavelet_norm))
end
%% Chi model
switch chi_model_type
    case 'empirical' 
        chi = (job_meta.s_rate/1e6)*(0:1:ns-1)'.*-2 + 19;

    case 'raw'
        
    case 'bootstrap'

    otherwise
        warning('Unexpected plot type. No plot created.');
end
%% Build operator matrix

% Build blanking matrix used to ensure the convolution operator matrix is correct
IGblank = spdiags(ones(ns,2*hns_wavelet+1),(-hns_wavelet:hns_wavelet),ns,ns);
IGblank = repmat(IGblank,1+job_meta.nvols,2);

% Tikhonov regularisation matrix % check tik order
%smooth = spdiags([-wsmooth*ones(2*ns,1) wsmooth*ones(2*ns,1)],[-1 0],2*ns,2*ns);

%smooth = spdiags([wsmooth*ones(2*ns,1) -2*wsmooth*ones(2*ns,1) wsmooth*ones(2*ns,1)],[-1 0 1],2*ns,2*ns);
smooth = spdiags([wsmooth*ones(2*ns,1) -2*wsmooth*ones(2*ns,1) wsmooth*ones(2*ns,1)],[-1 0 1],ns,ns);
smooth = [smooth, smooth];

% Extract the angle of the stack closest to the mean chi angle for this trace
% This is for the IG crossplot correlation constraint
theta = asind(sqrt(tand(mean(chi))));
[~, angle_index] = min(abs(input_angles-theta));

IGmatrix = build_operator(job_meta,input_angles,ns_wavelet,wavelet_norm,ns,hns_wavelet,angle_index,chi,eer_weight,IGblank,smooth);
%% Inversion loop

% Set first and last traces to loop over
first_iter = 1;
last_iter = ntraces;

% Begin inversion loop
tic
for kk = first_iter:last_iter    
    for ii = 1:job_meta.nvols
        data(:,ii) = traces{ii}(:,kk); % Read the angle stack data for the inversion of this trace
    end   
    data(isnan(data)) = 0;  % Set NaNs to zero
    fold = sum(data ~= 0,2); % Get angle fold
    data_tmp = data(:);  % Make temporary column vector from data
    
    data_zeros = abs(data_tmp) < abs(mean(data_tmp)/range(data_tmp)); % Find zones where data is zero (due to mute angle mute functions)
    % or close to zero
    data(reshape(data_zeros,[],job_meta.nvols)) = 0;
    data_zeros = logical([data_zeros;zeros(3*ns,1)]); 
    if plot_on == 1
        figure(6)
        subplot(1,2,1); spy(IGmatrix)
    end
    IGiter = IGmatrix;
    IGiter(data_zeros,:) = 0; % Set operator rows to zero if there are zeros in the data vector
    
    if plot_on == 1
        figure(6)
        subplot(1,2,2); spy(IGiter)
    end
    
    % Model intercept = [near], model gradient = [far-near]/[far_angle - near_angle]
    model_tmp = zeros(2,ns);
    for ii = 1:ns
        model_op = [ones(job_meta.nvols,1),sind(input_angles').*sind(input_angles')];
        model_zeros = data(ii,:) == 0;
        model_op(model_zeros,:) = 0;
        model_tmp(:,ii) = model_op\data(ii,:)';
    end
    model_tmp = model_tmp';
    %[traces{1}(:,1) traces{2}(:,1) traces{3}(:,1) model_tmp]; 
    if output_std == 1;
        Imodel(:,kk) = model_tmp(:,1)/norm(model_tmp(:,1)); %data(:,1)/norm(data(:,1));
        Gmodel(:,kk) = model_tmp(:,2)/norm(model_tmp(:,1)); %-Imodel./tand(chi);
        model = [Imodel(:,kk);Gmodel(:,kk)];
    else
        Imodel = model_tmp(:,1)/norm(model_tmp(:,1)); %data(:,1)/norm(data(:,1));
        Gmodel = model_tmp(:,2)/norm(model_tmp(:,1)); %-Imodel./tand(chi);
        %model = [Imodel;Gmodel];
    end
    % Set NaNs to zero
    model(isnan(model)) = 0;
    
    % Make the data a column vector and add zeros on the end for the EER constraint and the Tikhonov regularisation
    %data = [data(:);zeros(3*ns,1)];
    data = [data(:);zeros(2*ns,1)];
    % Do the inversion
    if background == 1
        [ava(:,kk),~] = lsqr(IGiter,data,tol,iter,[],[],model);   
        % try a different solving method?
    else     
        [ava(:,kk),~] = lsqr(IGiter,data,tol,iter,[],[]);
    end
    %[ava(:,kk),~] = lsqr(IGiter,data,1e-2,100,[],[],model);
    ava([fold;fold]==0,kk)=0;
    
    % Estimate the R^2 confidence in the result. This is the variance ratio:
    % 1-[(sum(data-Gm)^2)/(sum(data-mean)^2)], where the inversion was
    % solving data = Gm.
    data = reshape(data(1:ns*job_meta.nvols,:),[],job_meta.nvols);
    digi_confidence(:,kk) = 1-(sum((data-reshape(IGiter(1:job_meta.nvols*ns,:)*ava(:,kk),[],job_meta.nvols)).^2,2)./sum(bsxfun(@minus,data,sum(data,2)./fold).^2,2));
    
    % Give a status report, add this to a log file
    fprintf('Completed trace %d of %d\n',kk-first_iter+1,last_iter-first_iter+1)
end
toc

%% Save results
digi_intercept = [ava(1:ns,:);zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
digi_gradient = [ava(1+ns:end,:);zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
digi_minimum_energy_eer_projection = [bsxfun(@times,ava(1:ns,:),cosd(chi))+bsxfun(@times,ava(1+ns:end,:),sind(chi));zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];

digi_confidence(digi_confidence<0)=0;
digi_confidence(digi_confidence>1)=1;
digi_confidence(isnan(digi_confidence))=0;
digi_confidence = [digi_confidence;zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
if output_std == 1;
    std_intercept = [Imodel;zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
    std_gradient = [Gmodel;zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
    std_minimum_energy_eer_projection = [bsxfun(@times,Imodel,cosd(chi))+bsxfun(@times,Gmodel,sind(chi));zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
end

% Unflatten data
for kk = 1:length(wb_idx)
    digi_intercept(:,kk) = circshift(digi_intercept(:,kk),wb_idx(kk));
    digi_gradient(:,kk) = circshift(digi_gradient(:,kk),wb_idx(kk));
    digi_minimum_energy_eer_projection(:,kk) = circshift(digi_minimum_energy_eer_projection(:,kk),wb_idx(kk));
    if output_std == 1;
        std_intercept(:,kk) = circshift(std_intercept(:,kk),wb_idx(kk));
        std_gradient(:,kk) = circshift(std_gradient(:,kk),wb_idx(kk));
        std_minimum_energy_eer_projection(:,kk) = circshift(std_minimum_energy_eer_projection(:,kk),wb_idx(kk));
    end
end

% Save outputs into correct structure to be written to SEGY.
results_out{1,1} = 'Meta data for output files';
results_out{1,2}{1,1} = ilxl_read{vol_index_wb};
results_out{1,2}{2,1} = uint32(zeros(size(traces{1},2),1));
results_out{2,1} = 'digi_intercept';
results_out{2,2} = digi_intercept;
results_out{3,1} = 'digi_gradient';
results_out{3,2} = digi_gradient;
results_out{4,1} = 'digi_minimum_energy_eer_projection'; 
results_out{4,2} = digi_minimum_energy_eer_projection;
results_out{5,1} = 'digi_confidence'; 
results_out{5,2} = digi_confidence;
% output standard intercept and gradient result
if output_std == 1;
    results_out{6,1} = 'std_intercept';
    results_out{6,2} = std_intercept;
    results_out{7,1} = 'std_gradient';
    results_out{7,2} = std_gradient;
    results_out{8,1} = 'std_minimum_energy_eer_projection';
    results_out{8,2} = std_minimum_energy_eer_projection;
end

% check segy write functions - many different versions now!
if exist(strcat(job_meta.output_dir,'digi_results/','intercept_gradient_segy/'),'dir') == 0
    output_dir = strcat(job_meta.output_dir,'digi_results/','intercept_gradient_segy/');
    mkdir(output_dir);    
else
    output_dir = strcat(job_meta.output_dir,'digi_results/','intercept_gradient_segy/');
end

node_segy_write(results_out,i_block,job_meta.s_rate/1000,output_dir)

end

function [IGmatrix] = build_operator(job_meta,input_angles,ns_wavelet,wavelet_tmp,ns,hns_wavelet,angle_index,chi,alpha,IGblank,smooth)
    % Build operator for inversion
    for ii = 1:job_meta.nvols
        Iwavelet_interp(:,1+(ii-1)*ns_wavelet:ii*ns_wavelet) = wavelet_tmp{ii}';
        Gwavelet_interp(:,1+(ii-1)*ns_wavelet:ii*ns_wavelet) = wavelet_tmp{ii}'*(sind(input_angles(ii)).*sind(input_angles(ii)));
        %Iwavelet_interp(:,1+(ii-1)*ns_wavelet:ii*ns_wavelet) = interp1(wavelet_z_grid,wavelet_tmp{ii}',start_interp:1:end_interp,'linear','extrap');
        %Gwavelet_interp(:,1+(ii-1)*ns_wavelet:ii*ns_wavelet) = interp1(wavelet_z_grid,...
        %    wavelet_tmp{ii}'*(sind(input_angles(ii)).*sind(input_angles(ii))),start_interp:1:end_interp,'linear','extrap');
    end

    IGdiagnals = sort(reshape([(-hns_wavelet:hns_wavelet)',bsxfun(@plus,(-hns_wavelet:hns_wavelet)',(-ns:-ns:-ns*(job_meta.nvols-1)))],1,[]),'descend');

    Imatrix = spdiags(Iwavelet_interp,IGdiagnals,ns*job_meta.nvols,ns);
    Gmatrix = spdiags(Gwavelet_interp,IGdiagnals,ns*job_meta.nvols,ns);

    EERmatrix = alpha*[bsxfun(@times,Imatrix(1+ns*(angle_index-1):+ns*angle_index,:),cosd(chi)),bsxfun(@times,Imatrix(1+ns*(angle_index-1):+ns*angle_index,:),sind(chi))];

    IGmatrix = [[Imatrix,Gmatrix;EERmatrix].*IGblank;smooth];
end

function wavelet_interp = interpolate_wavelets(wavelet_z_grid,wavelet,start_interp,end_interp)
%     start_interp = min(wavelet_z_grid)-mode(diff(wavelet_z_grid'));
%     end_interp = max(wavelet_z_grid)+mode(diff(wavelet_z_grid'));
    wavelet_interp = interp1(wavelet_z_grid,wavelet',start_interp:1:end_interp,'linear','extrap');
end

function wavelet_norm = wavelet_rms_norm(job_meta,wavelet_interp,wavelet_z_grid,vol_index)
    % Normalise the wavelets to have constant energy w.r.t. angle. The energy
    % is set to that of the nearest angle wavelets. Wavelet energy still varies
    % w.r.t. time.
    A = cell2mat(wavelet_interp);
    %B = sqrt(sum(A.^2));
    B = sqrt(max(A.^2));
    C = reshape(B',length(wavelet_z_grid),[]);
    D = C(:,vol_index);
    for ii=1:job_meta.nvols
        E = A(:,1+(ii-1)*length(wavelet_z_grid):ii*length(wavelet_z_grid));
        F = bsxfun(@rdivide,bsxfun(@times,E,D'),sqrt(sum(E.^2)));
        wavelet_norm{ii} = F;
    end
end
