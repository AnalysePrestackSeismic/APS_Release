function anomalous_body_connector(job_meta_path,vol_index,anomalous_threshold,connectivity,area_thresh_lt,area_thresh_gt,join_block_id)
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
%% Function Defination: This operates on the binary files created by seismic anomaly spotter and joins up geobodies connecting similar anomalies.

% INPUTS:
%   job_meta_path:         Path to job meta file
%   vol_index:             Index of which volume to use
%   anomalous_threshold:   Make this number between 0 and 1 (hint: 0.9 means it will choose only the top 10 percentile points)
%   connectivity:          Prameter controlling connectivity6 (have common faces) , 18 (have at least one common edge ) or 26 (have at least one common vertex)
%   area_threshold_lt:     Threshold are to output horizons lesser than this value
%   area_threshold_gt:     Threshold are to output horizons greater than this value
%   join_block_id:         The id of the block that will combine results, set to 0 if no

% OUTPUTS:
%   void output

% WRITES TO DISK:
%       A binary anomaly volume 1 if anomaly greater than threshold, 0 if
%       anomaly lesser than threshold (file stored in binary format )

%       A series of horizons tops and base in a signature stamped horizon
%       directory. Only horizons greater than provided area_thresh are written out
%       ( as text files and images)

% RUNNING INSTRUCTION: Currently this is not submitted through slurm. So locally run on your
% machine preferably a powerful machine

%% Check if image processing tool box is available
write_vol=0 ; % Make this parameter 0 if you dont want write out the connectted bodies as volume 

[TF errmsg] = license('checkout','image_toolbox');                              %Check if emage processing toolbox is available
if TF == 0
    fprintf('\n No image_toolbox license available (James probably has it...) \n');
else
    make_blocky = 0;
    job_meta = load(job_meta_path);
    
    pkey_inc_mode = mode(job_meta.pkey_inc);
    skey_inc_mode = mode(job_meta.skey_inc);
    
    pkeyn = 1+((job_meta.pkey_max(str2double(vol_index))-job_meta.pkey_min(str2double(vol_index)))...
        /job_meta.pkey_inc(str2double(vol_index)));
    skeyn = 1+((job_meta.skey_max(str2double(vol_index))-job_meta.skey_min(str2double(vol_index)))...
        /job_meta.skey_inc(str2double(vol_index)));
    
    anomalous_threshold = str2num(anomalous_threshold);
    connectivity = str2num(connectivity);
    join_block_id = str2num(join_block_id);
    area_thresh_gt = str2num(area_thresh_gt);
    area_thresh_lt = str2num(area_thresh_lt);
    
    %% Load traces
    if isfield(job_meta, 'anom_result')
        %n_slices = job_meta.n_samples{1};
        %n_slices = 2330;
        
        % Load the water bottom file
        fprintf('Loading water bottom\n');
        wb_idx = dlmread(job_meta.wb_path); % import of ascii horizon
        wb_idx = sortrows(wb_idx,[1 2]); % sort inlines and crosslines because water bottom written out in block order
        
        wb_idx = wb_idx(wb_idx(:,1) > job_meta.pkey_min(1) &  wb_idx(:,1) < job_meta.pkey_max(1) & wb_idx(:,2) > job_meta.skey_min(1) & wb_idx(:,2) < job_meta.skey_max(1),:);
        
        n_iline = (wb_idx(:,1)-job_meta.pkey_min(1))/pkey_inc_mode+1;
        n_xline = (wb_idx(:,2)-job_meta.skey_min(1))/skey_inc_mode+1;
        
        lin_ind = ((n_iline-1).*skeyn)+n_xline;
        wb_ind = ones(1,pkeyn*skeyn);
        wb_ind(lin_ind) = floor(wb_idx(:,3)./(job_meta.s_rate/1000));
        wb_ind = medfilt1(wb_ind,3);
        fprintf('Loaded water bottom\n');
        clear wb_idx
        
        %     This section of code didn't work if the folder names had _ in them.
        %     file_parts = regexp(job_meta.anom_result, '\_', 'split');
        %     ind = strfind(file_parts{7},'samp');
        file_parts = regexp(job_meta.anom_result, 'sas_combine', 'split');
        file_parts = regexp(file_parts{2}, '\_', 'split');
        ind = strfind(file_parts{2},'samp');
        %n_slices = 1374;
        n_slices = str2double(file_parts{2}(ind+4:end));
        n_slices_unflat = n_slices+max(wb_ind);
        I = false(n_slices_unflat,pkeyn*skeyn); % binary volume
        block_sz = 50000; % read 50000 traces at a trace
        
        win_sub = ones(n_slices_unflat,block_sz);
        win_sub = cumsum(win_sub,1);
        
        fid_read = fopen(job_meta.anom_result,'r');
        for ii_block = 1:block_sz:pkeyn*skeyn
            maxblock = ii_block + block_sz - 1;
            
            if maxblock >  pkeyn*skeyn
                maxblock = pkeyn*skeyn;
                block_sz = (maxblock - ii_block+1);
                
                %win_sub = ones(n_slices_unflat,block_sz);
                %win_sub = cumsum(win_sub,1);
                win_sub = win_sub(:,1:block_sz);
            end
            traces_tmp = fread(fid_read,[n_slices,block_sz],'float32=>float32');
            %win_ind = true(n_slices_unflat,block_sz);
            %bsxfun(@plus,wb_ind(ii_block:maxblock),(0:n_slices_unflat-max(wb_ind(ii_block:maxblock)))');
            win_ind = true(n_slices_unflat,block_sz);
            win_ind(bsxfun(@lt,win_sub,wb_ind(ii_block:maxblock)) == 1) = 0;
            win_ind(bsxfun(@gt,win_sub,wb_ind(ii_block:maxblock)+n_slices-1) == 1) = 0;
            
            %win_sub_log = true(n_slices_unflat,block_sz);
            %win_sub_log(win_ind) = 0;
            %win_ind = bsxfun(@gt,win_sub,wb_ind(ii_block:maxblock)+n_slices-1);
            %win_sub_log(win_ind) = 0;
            %win_sub = win_sub > 0;
            traces_tmp_unflat = false(n_slices_unflat,block_sz);
            %traces_tmp = traces_tmp >= anomalous_threshold;
            traces_tmp_unflat(win_ind) = traces_tmp >= anomalous_threshold;
            I(:,ii_block:maxblock) = traces_tmp_unflat; %traces_tmp_unflat;
            %         if ii_block == 7350001
            %           fprintf('Comp');
            %         end
            fprintf('Completed reading traces %d to %d of %d\n',ii_block,maxblock,pkeyn*skeyn);
            clear traces_tmp
        end
        fclose(fid_read);
        
        I = I(1:n_slices,:);
        %     fid = fopen(strcat(job_meta.output_dir,sprintf('%s_unflat_anomalies_result_%s.bin',num2str(anomalous_threshold),'I')),'w');
        %     fwrite(fid,I,'float32');
        %     fclose(fid);
    else
        loopfin = size(job_meta.liveblocks,1);
        lpi = 1;
        while lpi <= loopfin
            i_block = job_meta.liveblocks(lpi);
            [~, traces, ilxl_read, ~] = ...
                node_segy_read(job_meta_path,vol_index,num2str(i_block));
            traces = [zeros(1,size(traces,2)); traces(2:end,:)];
            
            [n_slices,~] = size(traces);
            
            n_iline = (ilxl_read(:,1)-job_meta.pkey_min(str2double(vol_index)))/pkey_inc_mode+1;
            n_xline = (ilxl_read(:,2)-job_meta.skey_min(str2double(vol_index)))/skey_inc_mode+1;
            lin_ind = ((n_iline-1).*skeyn)+n_xline;
            
            if lpi == 1;
                I = false(n_slices,pkeyn*skeyn);
            end
            
            I(1:n_slices,double(lin_ind)) = traces >= anomalous_threshold;
            fprintf('-- Block %d --\n',lpi);
            lpi = lpi + 1;
            % fprintf('-- Block %d --\n',lpi);
            % Save water bottom pick
            
        end
    end
    %end
    %I = I(1:n_slices,:);
    if make_blocky == 1
        filtw = [1 2 3 2 1]/9;
        for i_trace = 1:1:size(I,2);
            I(:,i_trace) = conv(single(I(:,i_trace)),filtw,'same');
        end
        I = logical(I);
    end
    clear lin_ind n_iline n_xline traces_tmp_unflat wb_ind win_ind win_sub win_sub_log file_parts
    
    I = reshape(I,n_slices,skeyn,pkeyn);
    
    % run connectivity and stats algorithm
    % fprintf('-- Block %d --\n',process_positions.block_id);
    fprintf('Finding connected bodies\n');
    CC = bwconncomp(I,connectivity);
    fprintf('Completed and found %d\n',CC.NumObjects);
    STATS = regionprops(CC,'Area','BoundingBox','Centroid','PixelList','PixelIdxList');
    
    NumObjects = CC.NumObjects;
    if NumObjects <= intmax('uint8')
        dataType = 'uint8';
    elseif NumObjects <= intmax('uint16')
        dataType = 'uint16';
    elseif NumObjects <= intmax('uint32')
        dataType = 'uint32';
    else
        dataType = 'double';
    end    
    
    A = cell2mat(CC.PixelIdxList');
    
    clear CC
    clear I;
    n_vol = 40; % need work this out dynamically
    length_I = n_slices*skeyn*pkeyn;
    for ii_n_vol = 1:n_vol-1;
        length_I_vol(ii_n_vol) = floor(length_I/n_vol);
    end
    length_I_vol(n_vol) = length_I - (n_vol-1)*floor(length_I/n_vol);
    
    ind = 1; ii = 1; while ii < size(A,1); A(ii:ii+STATS(ind).Area-1,2) = STATS(ind).Area; ii = ii+STATS(ind).Area; ind = ind+1; end
    
    start_il = job_meta.pkey_min;
    start_xl = job_meta.skey_min;
    n_xls = skeyn;
    n_samps = n_slices;
    
    %fid = fopen(strcat(job_meta.output_dir,sprintf('%s_anomalies_result_%s_sil%s_sxl%s_nxl%s_samp%s.bin',num2str(anomalous_threshold),'I_id',num2str(start_il),num2str(start_xl),num2str(n_xls),num2str(n_samps))),'a');
    A = sortrows(A,1);
    
    %last_ind = 0;
    %--------This is for writing the volume of anomalous bodies, comment
    %out if you dont need them-----------------
    if write_vol==1
        id = 1;
        fprintf('Writing out volume...\n');
        for ii_vol = 1:n_vol
            fid = fopen(strcat(job_meta.output_dir,sprintf('%s_anomalies_result_%s_sil%s_sxl%s_nxl%s_samp%s.bin',num2str(anomalous_threshold),'I_vol',num2str(start_il),num2str(start_xl),num2str(n_xls),num2str(n_samps))),'a');
            I_vol = zeros(length_I_vol(ii_vol),1,dataType);
            %cur_ind = length_I_vol(ii_vol)+last_ind;
            ind = A(:,1) >= ((ii_vol-1)*length_I_vol(ii_vol))+1 & A(:,1) <= ii_vol*length_I_vol(ii_vol);
            ind_lin = A(ind,1) - (ii_vol-1)*length_I_vol(ii_vol);
            I_vol(ind_lin) = A(ind,2);
            fwrite(fid,I_vol,'float32');
            fclose(fid);
            %last_ind = last_ind + length_I_vol(ii_vol);
            
            fprintf('Completed %d of %d\n',ii_vol,n_vol);
            
            %         fid = fopen(strcat(job_meta.output_dir,sprintf('%s_anomalies_result_%s_sil%s_sxl%s_nxl%s_samp%s.bin',num2str(anomalous_threshold),'I_id',num2str(start_il),num2str(start_xl),num2str(n_xls),num2str(n_samps))),'a');
            %         I_vol = zeros(length_I_vol(ii_vol),n_vol,dataType);
            %         ind = A(:,1) >= ((ii_vol-1)*length_I_vol(ii_vol))+1 & A(:,1) <= ii_vol*length_I_vol(ii_vol);
            %         ind_lin = A(ind,1) - (ii_vol-1)*length_I_vol(ii_vol);
            %         I_vol(ind_lin) = A(ind,2);
            %         fwrite(fid,I_vol,'float32');
            %         fclose(fid);
            
            id = id + 1;
        end
    end
    
    fprintf('Creating surfaces...\n'); % could add pause and it wait for arguments
    % Create surfaces
    str_date = date;
    str_date = regexprep(str_date, '-', '');
    out_dir_top = strcat(job_meta.output_dir,'/horizons_top',str_date,'_',num2str(anomalous_threshold),'_',num2str(area_thresh_lt),'_',num2str(area_thresh_gt));
    out_dir_bot = strcat(job_meta.output_dir,'/horizons_bot',str_date,'_',num2str(anomalous_threshold),'_',num2str(area_thresh_lt),'_',num2str(area_thresh_gt));
    out_dir_img = strcat(job_meta.output_dir,'/horizons_img',str_date,'_',num2str(anomalous_threshold),'_',num2str(area_thresh_lt),'_',num2str(area_thresh_gt));
    out_dir_pol = strcat(job_meta.output_dir,'/horizons_poly',str_date,'_',num2str(anomalous_threshold),'_',num2str(area_thresh_lt),'_',num2str(area_thresh_gt));
    mkdir(out_dir_top);
    mkdir(out_dir_bot);
    mkdir(out_dir_img);
    mkdir(out_dir_pol);
    
    pkeys = job_meta.pkey_min:job_meta.pkey_inc:job_meta.pkey_max;
    skeys = job_meta.skey_min:job_meta.skey_inc:job_meta.skey_max;
    
    for ii = 1:1:NumObjects
        if STATS(ii).Area > area_thresh_gt && STATS(ii).Area < area_thresh_lt
            fprintf('Creating surfaces for body %d of %d\n',ii,NumObjects);
            z_col = 2;
            
            skey_col = 1;
            pkey_col = 3;
            
            pos = STATS(ii).PixelList;
            sort_z_as = sortrows(pos,z_col);
            sort_z_ds = sortrows(pos,-z_col);
            pos_uniq = unique(pos(:,[pkey_col skey_col]),'rows');
            % col 1 skey
            % col 2 z
            % col 3 pekyn
            pos_body_surface_top(:,1) = pkeys(pos_uniq(:,1));
            pos_body_surface_top(:,2) = skeys(pos_uniq(:,2));
            
            [~,Locb] = ismember(pos_uniq,sort_z_ds(:,[pkey_col skey_col]),'rows');
            pos_body_surface_top(:,3) = (sort_z_ds(Locb,2)-50).*job_meta.s_rate/1000;
            
            pos_body_surface_bot(:,1) = pkeys(pos_uniq(:,1));
            pos_body_surface_bot(:,2) = skeys(pos_uniq(:,2));
            [~,Locb] = ismember(pos_uniq,sort_z_as(:,[3 1]),'rows');
            pos_body_surface_bot(:,3) = (sort_z_as(Locb,2)-50).*job_meta.s_rate/1000;
            z_mode = mode(pos_body_surface_top(:,3));
            % add mean depth to filename
            top_asc = strcat(out_dir_top,'/anom_surface_body_z',num2str(z_mode),'_area',num2str(STATS(ii).Area),'_id',num2str(ii),'_top.txt');
            bot_asc = strcat(out_dir_bot,'/anom_surface_body_z',num2str(z_mode),'_area',num2str(STATS(ii).Area),'_id',num2str(ii),'_bot.txt');
            
            area = repmat(STATS(ii).Area,size(pos_body_surface_bot(:,3),1),1);
            id = repmat(ii,size(pos_body_surface_bot(:,3),1),1);
            
            dlmwrite(top_asc,[pos_body_surface_bot area id],'delimiter','\t','precision', 8)
            dlmwrite(bot_asc,[pos_body_surface_top area id],'delimiter','\t','precision', 8)
            
            min_il = min(pos_body_surface_top(:,1));
            max_il = max(pos_body_surface_top(:,1));
            min_xl = min(pos_body_surface_top(:,2));
            max_xl = max(pos_body_surface_top(:,2));            
            
            iii = 1;
            for il_ii = min_il:job_meta.pkey_inc:max_il
                poly_min(iii,:) = [il_ii,min(pos_body_surface_top((pos_body_surface_top(:,1) == il_ii),2))];
                poly_max(iii,:) = [il_ii,max(pos_body_surface_top((pos_body_surface_top(:,1) == il_ii),2))];
                iii = iii + 1;
            end
            
            for xl_ii = min_xl:job_meta.skey_inc:max_xl
                poly_min(iii,:) = [min(pos_body_surface_top((pos_body_surface_top(:,2) == xl_ii),1)),xl_ii];
                poly_max(iii,:) = [max(pos_body_surface_top((pos_body_surface_top(:,2) == xl_ii),1)),xl_ii];
                iii = iii + 1;
            end
            
            poly = [poly_min;poly_max];
            poly = unique(poly,'rows');
            %poly = sortrows(poly);
            
            %poly_area = 
            %mid_il = min_il+floor((max_il-min_il)/2);
            %mid_xl = min_xl+floor((max_xl-min_xl)/2);
                        
            %[theta,~] = cart2pol(poly(:,1)-mid_il,poly(:,2)-mid_xl);
            %poly(:,3) = theta;
            %poly = sortrows(poly,3);
            poly(:,3) = pos_body_surface_top(ismember(pos_body_surface_top(:,1:2),poly(:,1:2),'rows'),3);
            pol_asc = strcat(out_dir_pol,'/anom_poly_z',num2str(z_mode),'_area',num2str(STATS(ii).Area),'_id',num2str(ii),'.txt');
            dlmwrite(pol_asc,poly,'delimiter','\t','precision', 8)         
           
            
            clear poly* pos sort_z_as sort_z_ds pos_uniq Locb pos_body_surface_top pos_body_surface_bot
            
        end
    end
    
end
end
