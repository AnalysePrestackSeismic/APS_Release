function [] = int_grad_inv_proj_lite_js(seismic_mat_path,i_block,n_blocks,wavelet_mat_path,output_dir)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% INTERCEPT GRADIENT INVERSION
% Date: 15 October 2012
% Authors: Jonathan Edgar and James Selvage

    i_block = str2double(i_block);
    n_blocks = str2double(n_blocks);

    n_blocks = 1;
    tmp_seismic = fread(fopen(seismic_mat_path,'r'),'double');
    
    seismic.n_samples = tmp_seismic(1,1);
    seismic.trace_ilxl_bytes = reshape(tmp_seismic(2:end),3,[])';
    seismic.n_traces = size(seismic.trace_ilxl_bytes,1);
    
    seismic.srate = 0.004;
    seismic.nfiles = 6;
    seismic.angles = [7.5,12.5,17.5,22.5,27.5,32.5];
    
    load(wavelet_mat_path);
    
    wavelet_z_grid = unique(wavelet{1}(1,:)); % gets sample of midpoint of each window
    [~,wavelet_tmp_idx] = ismember(wavelet_z_grid,fliplr(wavelet{1}(1,:)));
    wavelet_tmp_idx = size(wavelet{1},2)-wavelet_tmp_idx+1;
    
    dec = 0;
    output_std = 0;

    for ii = 1:seismic.nfiles
        if ii == 1
            seismic.filepath = '/segy/KEN/2012_L10ab_final/filtered_angle_stacks/113j02_filtered_angle_stack_05-10_SEGY_disk.segy';
%             seismic.filepath = '/data/TZA/dtect/TZA_Kusini_Outboard/matlab/kusini_2012_113j05_05-10_angle_stk_merge_il4634';
            [traces{ii}, ilxl_read{ii}] = node_segy_read_traces_lite(seismic,i_block,n_blocks,dec);
        else
            seismic.filepath = sprintf('/segy/KEN/2012_L10ab_final/filtered_angle_stacks/113j02_filtered_angle_stack_%d-%d_SEGY_disk.segy',ii*5,(ii+1)*5);
%             seismic.filepath = sprintf('/data/TZA/dtect/TZA_Kusini_Outboard/matlab/kusini_2012_113j05_%d-%d_angle_stk_merge_il4634',ii*5,(ii+1)*5);
            [traces{ii}, ilxl_read{ii}] = node_segy_read_traces_lite(seismic,i_block,n_blocks,dec);
        end
    end 
    
    wb_idx = water_bottom_flatten_lite(traces{4}); %Runs water bottom flatten on mid offsets?
    
    for ii = 1:seismic.nfiles
        for kk = 1:length(wb_idx)
            traces{ii}(:,kk) = circshift(traces{ii}(:,kk),-wb_idx(kk));
            traces{ii}(end-wb_idx(kk):end,kk) = 0;
        end
        traces{ii} = traces{ii}(1:max(wavelet_z_grid)+min(wavelet_z_grid),:);
    end

    % Initial some variables
    ns_wavelet = size(wavelet{1},1)-1;
    hns_wavelet = floor(ns_wavelet/2);
    [ns ntraces] = size(traces{1});
    alpha = 1; % Weight for EER constraint
    
    % Build blanking matrix used to ensure the convolution operator matrix is correct
    IGblank = spdiags(ones(ns,2*hns_wavelet+1),(-hns_wavelet:hns_wavelet),ns,ns);
    IGblank = repmat(IGblank,1+seismic.nfiles,2);
    
    % Tikhonov regularisation weight
    wsmooth = 0;
    % Tikhonov regularisation matrix
    smooth = spdiags([-wsmooth*ones(2*ns,1) 2*wsmooth*ones(2*ns,1) -wsmooth*ones(2*ns,1)],[-1 0 1],2*ns,2*ns);
    
    % Build the seabed taper
%     taper_length = 50;
%     taper = (1+cos(pi:-pi/(taper_length-1):0))/2;
%     taper = [taper(:);ones(ns-taper_length,1)];

%% Inversion loop
    
    % Set first and last traces to loop over
    first_iter = 1;
    last_iter = ntraces;
    
    % Begin inversion loop
    tic
    for kk = first_iter:last_iter

        % Read the angle stack data for the inversion of this trace
        for ii = 1:seismic.nfiles
            data(:,ii) = traces{ii}(:,kk);     
        end
        
%         data = bsxfun(@times,data,taper);
        
        % Build the minimum energy chi angle from the chi model for this trace location
%        chi(:,kk) = (1:1:ns)'.*chi_model_grid(chi_finder(:,kk),4) + chi_model_grid(chi_finder(:,kk),3);
        chi = seismic.srate*(0:1:ns-1)'.*-1 + 20;
        
        % Extract the angle of the stack closest to the mean chi angle for this trace
        theta = asind(sqrt(tand(mean(chi))));
        [~, angle_index] = min(abs(seismic.angles-theta));
        
        if kk == first_iter
            for ii=1:seismic.nfiles
                wavelet_tmp{ii} = wavelet{ii}(2:end,:);
            end
            [IGmatrix] = build_operator(seismic,ns_wavelet,wavelet_z_grid,wavelet_tmp,wavelet_tmp_idx,ns,hns_wavelet,angle_index,chi,alpha,IGblank,smooth);
        end
                     
        % Set NaNs to zero
        data(isnan(data)) = 0;
        
        % Make temporary column vector from data
        data_tmp = data(:);
        
        % Find zones where data is zero (due to mute angle mute functions)
        data_zeros = data_tmp == 0;
        data_zeros = logical([data_zeros;zeros(3*ns,1)]);
        
        % Set operator rows to zero if there are zeros in the data vector
        IGmatrix(data_zeros,:) = 0;
            
        % Model intercept = [near], model gradient = [far-near]/[far_angle - near_angle]
        model_tmp = zeros(2,ns);
        for ii = 1:ns
            model_op = [ones(seismic.nfiles,1),sind(seismic.angles').*sind(seismic.angles')];
            model_zeros = data(ii,:) == 0;
            model_op(model_zeros,:) = 0;
            model_tmp(:,ii) = model_op\data(ii,:)';
        end
        model_tmp = model_tmp';
        if output_std == 1;
            Imodel(:,kk) = model_tmp(:,1)/norm(model_tmp(:,1)); %data(:,1)/norm(data(:,1));
            Gmodel(:,kk) = model_tmp(:,2)/norm(model_tmp(:,1)); %-Imodel./tand(chi);
            Gmodel_test(:,kk) = model_tmp(:,2)/norm(model_tmp(:,2)); %-Imodel./tand(chi);
            model = [Imodel(:,kk);Gmodel(:,kk)];
        else
            Imodel = model_tmp(:,1)/norm(model_tmp(:,1)); %data(:,1)/norm(data(:,1));
            Gmodel = model_tmp(:,2)/norm(model_tmp(:,1)); %-Imodel./tand(chi);
            Gmodel_test(:,kk) = model_tmp(:,2)/norm(model_tmp(:,2)); %-Imodel./tand(chi);
            model = [Imodel;Gmodel];
        end
        % Set NaNs to zero
        model(isnan(model)) = 0;
        
        % Make the data a column vector and add zeros on the end for the EER constraint and the Tikhonov regularisation
        data = [data(:);zeros(3*ns,1)];

        % Do the inversion
        [ava(:,kk),~] = lsqr(IGmatrix,data,1e-2,100,[],[],model);

        % Clear the data ready for the next trace
        clearvars data
        
        % Give a status report
        fprintf('Completed trace %d of %d\n',kk-first_iter+1,last_iter-first_iter+1)
    end
    toc
    
    % Make outputs to be saved to .mat files
    digi_intercept = [ava(1:ns,:);zeros(seismic.n_samples-ns,ntraces)];
    digi_gradient = [ava(1+ns:end,:);zeros(seismic.n_samples-ns,ntraces)];
    digi_minimum_energy_eer_projection = [bsxfun(@times,ava(1:ns,:),cosd(chi))+bsxfun(@times,ava(1+ns:end,:),sind(chi));zeros(seismic.n_samples-ns,ntraces)];
    if output_std == 1;
        std_intercept = [Imodel;zeros(seismic.n_samples-ns,ntraces)];
        std_gradient = [Gmodel;zeros(seismic.n_samples-ns,ntraces)];
        std_minimum_energy_eer_projection = [bsxfun(@times,Imodel,cosd(chi))+bsxfun(@times,Gmodel,sind(chi));zeros(seismic.n_samples-ns,ntraces)];
    end
    
    
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

    % Save outputs as .mat files
    results_out{1,1} = 'digi_intercept';
    results_out{1,2} = digi_intercept;
    results_out{2,1} = 'digi_gradient';
    results_out{2,2} = digi_gradient;
    results_out{3,1} = 'digi_minimum_energy_eer_projection';
    results_out{3,2} = digi_minimum_energy_eer_projection;
    results_out{4,1} = 'ilxl numbers';
    results_out{4,2} = ilxl_read;
    if output_std == 1;
        results_out{5,1} = 'std_intercept';
        results_out{5,2} = std_intercept;
        results_out{6,1} = 'std_gradient';
        results_out{6,2} = std_gradient;
        results_out{7,1} = 'std_minimum_energy_eer_projection';
        results_out{7,2} = std_minimum_energy_eer_projection;
    end
    
    save(strcat(output_dir,sprintf('%s_results_block_%d','digi',i_block)),'results_out','-v7.3');
end

%%
function [IGmatrix] = build_operator(seismic,ns_wavelet,wavelet_z_grid,wavelet_tmp,wavelet_tmp_idx,ns,hns_wavelet,angle_index,chi,alpha,IGblank,smooth)

    for ii = 1:seismic.nfiles
        Iwavelet_interp(:,1+(ii-1)*ns_wavelet:ii*ns_wavelet) = interp1(wavelet_z_grid,wavelet_tmp{ii}(:,wavelet_tmp_idx)',1:1:ns,'linear','extrap');
        Gwavelet_interp(:,1+(ii-1)*ns_wavelet:ii*ns_wavelet) = interp1(wavelet_z_grid,...
            wavelet_tmp{ii}(:,wavelet_tmp_idx)'*(sind(seismic.angles(ii)).*sind(seismic.angles(ii))),1:1:ns,'linear','extrap');
    end

    IGdiagnals = sort(reshape([(-hns_wavelet:hns_wavelet)',bsxfun(@plus,(-hns_wavelet:hns_wavelet)',(-ns:-ns:-ns*(seismic.nfiles-1)))],1,[]),'descend');

    Imatrix = spdiags(Iwavelet_interp,IGdiagnals,ns*seismic.nfiles,ns);
    Gmatrix = spdiags(Gwavelet_interp,IGdiagnals,ns*seismic.nfiles,ns);
    
    EERmatrix = alpha*[bsxfun(@times,Imatrix(1+ns*(angle_index-1):+ns*angle_index,:),cosd(chi)),bsxfun(@times,Imatrix(1+ns*(angle_index-1):+ns*angle_index,:),sind(chi))];

    IGmatrix = [[Imatrix,Gmatrix;EERmatrix].*IGblank;smooth];
end
