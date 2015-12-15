function int_grad_inv_proj(block_mat_all,block_mat, process_files_mat, wavelet_mat_dir,chi_mat_dir,output_dir)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% INTERCEPT GRADIENT INVERSION
% Date: 15 October 2012
% Authors: Jonathan Edgar and James Selvage
%

%% Load required mat files

    fprintf('Loading mat files\n') 
    % Load the seismic file names and locations required by the node
    process_files = load(process_files_mat);
    % Load all wavelets estimated from all angle stacks
    wavelet = load_mat_files(wavelet_mat_dir,'wavelet_block','wavelet');
    % Load the minimum energy chi model
    chi_model = load_mat_files(chi_mat_dir,'chi_model_block','chi_model');
    
%% Read the seismic data

    % Read data from all angle stacks into cell array
    fprintf('Loading seismic\n') 
    [traces, process_positions] = segy_read_multiple_files(block_mat_all,process_files_mat, block_mat);
    
%% Grid the chi model and work out the trace positions of the chi model and wavelets nearest to each inlne of the processing grid for the node
    
    % Grid the chi model to all traces of the processing grid for the node
    fprintf('Gridding chi model\n') 
    chi_model_grid = grid_chi_model(chi_model,process_files,process_positions);
    
    % Calculate the number of traces in the processing grid
    fprintf('Locating wavelet positions\n') 
    ntraces_grid = length(process_positions.ilxl_grid);
    
    wavelet_ils = unique(wavelet{1}(1,:))';
    wavelet_xls = unique(wavelet{1}(2,:))';
    
    wavelet_pos = sortrows([repmat(wavelet_ils,size(wavelet_xls,1),1),repmat(wavelet_xls,size(wavelet_ils,1),1)],[1,2]);

%     fprintf('Building chi model lookup matrix\n') 
%     % Work out the chi model column index of first crossline on the inline of the wavelet nearest to each inline of the processing grid
%      [~,nearest_chi_il_index] = min(abs(bsxfun(@minus,process_positions.ilxl_grid(:,1)',chi_model_grid(:,1))));
%      
%     for kk = 1:ntraces_grid
%         [~,nearest_chi_xl_index(1,kk)] = min(abs(process_positions.ilxl_grid(kk,2)-chi_model_grid(nearest_chi_il_index(kk):end,2)));
%         nearest_chi_ilxl_index(1,kk) = nearest_chi_xl_index(1,kk)+nearest_chi_il_index(1,kk)-1;
%         chi_finder(:,kk) = ismember(chi_model_grid(:,1:2),chi_model_grid(nearest_chi_ilxl_index(kk),1:2),'rows');
%         % Set up the wavelet z position grid   
%     end
    
    wavelet_z_grid = unique(wavelet{1}(3,:));
    [~,wavelet_tmp_idx] = ismember(wavelet_z_grid,fliplr(wavelet{1}(3,:)));
    wavelet_tmp_idx = size(wavelet{1},2)-wavelet_tmp_idx+1;
    
%% Initialise some variables, build the blanking matrix and Tikhonov regularisation
    
    % Initial some variables
    ns_wavelet = size(wavelet{1},1)-3;
    hns_wavelet = floor(ns_wavelet/2);
    [ns ntraces] = size(traces{1});
    alpha = 1; % Weight for EER constraint
    
    % Build blanking matrix used to ensure the convolution operator matrix is correct
    IGblank = spdiags(ones(ns,2*hns_wavelet+1),(-hns_wavelet:hns_wavelet),ns,ns);
    IGblank = repmat(IGblank,1+process_files.nfiles,2);
    
    % Tikhonov regularisation weight
    wsmooth = 0;
    % Tikhonov regularisation matrix
    smooth = spdiags([-wsmooth*ones(2*ns,1) 2*wsmooth*ones(2*ns,1) -wsmooth*ones(2*ns,1)],[-1 0 1],2*ns,2*ns);
    
    % Build the seabed taper
%     taper_length = 50;
%     taper = (1+cos(pi:-pi/(taper_length-1):0))/2;
%     taper = [taper(:);ones(ns-taper_length,1)];
    
    ilxl_wavelet_blend_distance = 25;
    ilxl_wavelet_blend_flag = 0;

%% Inversion loop
    
    % Set first and last traces to loop over
    first_iter = 1;
    last_iter = ntraces;
    
    % Begin inversion loop
    tic
    for kk = first_iter:last_iter

        % Read the angle stack data for the inversion of this trace
        for ii = 1:process_files.nfiles
            data(:,ii) = traces{ii}(:,kk);
        end
        
%         data = bsxfun(@times,data,taper);
        
        % Build the minimum energy chi angle from the chi model for this trace location
%        chi(:,kk) = (1:1:ns)'.*chi_model_grid(chi_finder(:,kk),4) + chi_model_grid(chi_finder(:,kk),3);
        chi(:,kk) = (1:1:ns)'.*chi_model_grid(kk,4) + chi_model_grid(kk,3);
        
        % Extract the angle of the stack closest to the mean chi angle for this trace
        theta = asind(sqrt(tand(mean(chi(:,kk)))));
        [~, angle_index] = min(abs(process_files.angle-theta));
  
        trace_pos = process_positions.ilxl_grid(kk,1:2);
        wavelet_distances = ((((wavelet_pos(:,1)-trace_pos(1,1))./process_positions.il_inc).^2)+(((wavelet_pos(:,2)-trace_pos(1,2))./process_positions.xl_inc).^2)).^0.5;
        [min_wavelet_distance,min_wavelet_distance_index] = min(wavelet_distances);
        min_wavelet_distance_index_tmp = zeros(size(wavelet_distances));
        min_wavelet_distance_index_tmp(min_wavelet_distance_index,1) = 1;
        min_wavelet_distance_index = logical(abs(min_wavelet_distance_index_tmp-1));
        
        if min(wavelet_distances(min_wavelet_distance_index) - min_wavelet_distance) <= ilxl_wavelet_blend_distance
            fprintf('Rebuilding operator\n')
            ilxl_wavelet_blend_flag = 1;
            [~,next_min_wavelet_distance_index] = min(wavelet_distances(min_wavelet_distance_index) - min_wavelet_distance);
            wavelet_pos_tmp = wavelet_pos(min_wavelet_distance_index,:);
            wavelet_distances_tmp = wavelet_distances(min_wavelet_distance_index,:);
            for ii=1:process_files.nfiles
                nearest_wavelet_tmp_index{ii} = ismember(wavelet{ii}(1:2,:)',wavelet_pos(~min_wavelet_distance_index,:),'rows');
                next_nearest_wavelet_tmp_index{ii} = ismember(wavelet{ii}(1:2,:)',wavelet_pos_tmp(next_min_wavelet_distance_index,:),'rows');
                ilxl_wavelet_blend_weight = 1-(wavelet_distances(~min_wavelet_distance_index,:)/((wavelet_distances(~min_wavelet_distance_index,:)+wavelet_distances_tmp(next_min_wavelet_distance_index,:))));
                wavelet_tmp{ii} = ilxl_wavelet_blend_weight.*wavelet{ii}(4:end,nearest_wavelet_tmp_index{ii}') + (1-ilxl_wavelet_blend_weight).*wavelet{ii}(4:end,next_nearest_wavelet_tmp_index{ii}');
            end
            [IGmatrix] = build_operator(process_files,ns_wavelet,wavelet_z_grid,wavelet_tmp,wavelet_tmp_idx,ns,hns_wavelet,angle_index,chi(:,kk),alpha,IGblank,smooth);
        elseif or(kk == first_iter,ilxl_wavelet_blend_flag == 1)
            fprintf('Rebuilding operator\n')
            for ii=1:process_files.nfiles
                nearest_wavelet_tmp_index{ii} = ismember(wavelet{ii}(1:2,:)',wavelet_pos(~min_wavelet_distance_index,:),'rows');
                wavelet_tmp{ii} = wavelet{ii}(4:end,nearest_wavelet_tmp_index{ii}');
            end
            [IGmatrix] = build_operator(process_files,ns_wavelet,wavelet_z_grid,wavelet_tmp,wavelet_tmp_idx,ns,hns_wavelet,angle_index,chi(:,kk),alpha,IGblank,smooth);
            ilxl_wavelet_blend_flag = 0;
        end
        
        % Set NaNs to zero
        data(isnan(data)) = 0;
            
        % Model intercept = [near], model gradient = [far-near]/[far_angle - near_angle]
        Imodel = data(:,1)/norm(data(:,1));
        Gmodel = -Imodel./tand(chi(:,kk));
        model = [Imodel;Gmodel];
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
    intercept = ava(1:ns,:);
    gradient = ava(1+ns:end,:);
    minimum_energy_eer_projection = intercept.*cosd(chi)+gradient.*sind(chi);

    % Save outputs as .mat files
    results_out{1,1} = 'intercept';
    results_out{1,2} = intercept;
    results_out{2,1} = 'gradient';
    results_out{2,2} = gradient;
    results_out{3,1} = 'chi';
    results_out{3,2} = chi;
    results_out{4,1} = 'minimum_energy_eer_projection';
    results_out{4,2} = minimum_energy_eer_projection;
    
    save(strcat(output_dir,sprintf('%s_results_block_%d',process_files.func_name,process_positions.block_id)),'results_out','-v7.3');
end

%%
function [IGmatrix] = build_operator(process_files,ns_wavelet,wavelet_z_grid,wavelet_tmp,wavelet_tmp_idx,ns,hns_wavelet,angle_index,chi,alpha,IGblank,smooth)

    for ii = 1:process_files.nfiles
        Iwavelet_interp(:,1+(ii-1)*ns_wavelet:ii*ns_wavelet) = interp1(wavelet_z_grid,wavelet_tmp{ii}(:,wavelet_tmp_idx)',1:1:ns,'linear','extrap');
        Gwavelet_interp(:,1+(ii-1)*ns_wavelet:ii*ns_wavelet) = interp1(wavelet_z_grid,...
            wavelet_tmp{ii}(:,wavelet_tmp_idx)'*(sind(process_files.angle(ii)).*sind(process_files.angle(ii))),1:1:ns,'linear','extrap');
    end

    IGdiagnals = sort(reshape([(-hns_wavelet:hns_wavelet)',bsxfun(@plus,(-hns_wavelet:hns_wavelet)',(-ns:-ns:-ns*(process_files.nfiles-1)))],1,[]),'descend');

    Imatrix = spdiags(Iwavelet_interp,IGdiagnals,ns*process_files.nfiles,ns);
    Gmatrix = spdiags(Gwavelet_interp,IGdiagnals,ns*process_files.nfiles,ns);
    
    EERmatrix = alpha*[bsxfun(@times,Imatrix(1+ns*(angle_index-1):+ns*angle_index,:),cosd(chi)),bsxfun(@times,Imatrix(1+ns*(angle_index-1):+ns*angle_index,:),sind(chi))];

    IGmatrix = [[Imatrix,Gmatrix;EERmatrix].*IGblank;smooth];
end

%%
function [load_out] = load_mat_files(input_dir,search_string,variable_name)

    start_location = pwd;
    cd(input_dir);

    % Figure out the number of files in the current directory
    [~,nfiles] = (system(sprintf('ls -B | grep %s | wc -l',search_string))); 
    nfiles=str2double(nfiles);
    
    for qq = 1:nfiles
        if qq == 1
            load(sprintf('%s_%d',search_string,qq),variable_name);
            load_tmp = eval(variable_name);
            if iscell(load_tmp)
                ncells = size(load_tmp,2);
                [nrows,ncols] = size(load_tmp{1});
                cell_flag = 1;
                load_tmp = [load_tmp{:}];
            else
                cell_flag = 0;
            end
                
        else
            load(sprintf('%s_%d',search_string,qq),variable_name);
            if iscell(eval(variable_name))
                load_tmp = [load_tmp;eval(strcat('[',variable_name,'{:}]'))];
            else
                load_tmp = [load_tmp;eval(variable_name)];
            end
        end
    end
    
    if cell_flag == 1
        for ii = 1:ncells
            load_out_tmp{ii} = load_tmp(:,1+(ii-1)*ncols:ii*ncols);
            for jj = 1:nfiles
                load_out{ii}(:,1+(jj-1)*ncols:jj*ncols) = load_out_tmp{ii}(1+(jj-1)*nrows:jj*nrows,:);
            end
        end
    else
        load_out = load_tmp;
    end

    cd(start_location);
end

%%
function [chi_model_grid] = grid_chi_model(chi_model,process_files,process_positions)

    [X,Y] = meshgrid(unique(chi_model(:,1)),unique((chi_model(:,2))));
    [Xi,Yi] = meshgrid(unique(process_files.processing_grid.ilxl_grid(:,1)),unique((process_files.processing_grid.ilxl_grid(:,2))));
    
    chi_model_intercept = reshape(chi_model(:,3),size(X,1),size(Y,2));
    chi_model_intercept_grid = interp2(X,Y,chi_model_intercept,Xi,Yi,'spline');
    
    chi_model_gradient = reshape(chi_model(:,4),size(X,1),size(Y,2));
    chi_model_gradient_grid = interp2(X,Y,chi_model_gradient,Xi,Yi,'spline');
    
    chi_model_grid = [process_files.processing_grid.ilxl_grid(:,1:2),chi_model_intercept_grid(:),chi_model_gradient_grid(:)];
    
    chi_model_il_grid_index = and(chi_model_grid(:,1) >= min(process_positions.ilxl_grid(:,1)),chi_model_grid(:,1) <= max(process_positions.ilxl_grid(:,1)));
    chi_model_grid = chi_model_grid(chi_model_il_grid_index,:);   
end