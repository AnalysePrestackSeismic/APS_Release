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
%   job_meta_path:      Path to job meta file
%   vol_index:          Index of which volume to use 
%   block_mat:          information about the traces to process in this block
%   process_files_mat:  structures for the files to be read in
%   anomalous_threshold:between 0 and 1
%   connectivity:       6 (have common faces) , 18 (have at least one common edge ) or 26 (have at least one common vertex)
%   area_threshold:     Throshold are to output horizons greter than this value
%   join_block_id - the id of the block that will combine results, set to 0 if no

% OUTPUTS:
%   void output
%   Writes to Disk:
%       A binary anomaly volume 1 if anomaly greater than threshold, 0 if
%       anomaly lesser than threshold (file stored in binary format )
%       A series of horizons tops and base in a signature stamped horizon
%       directory. Only horizons greater than provided area_thresh are written out
%       ( as text files and images)

%% Check if image processing tool box is available

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
        fprintf('Completed traces %d to %d of %d\n',ii_block,maxblock,pkeyn*skeyn);
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
%I = I(1:job_meta.n_samples{1},:,:);
%ns = size(I,1);
% for ii = 1:1:n_slices
%     I(ii,:,:) = reshape(anom_traces_logic(ii,:),nxl,nil);
% end

% run connectivity and stats algorithm
% fprintf('-- Block %d --\n',process_positions.block_id);
fprintf('Finding connected bodies\n');
CC = bwconncomp(I,connectivity);
fprintf('Completed and found %d\n',CC.NumObjects);
STATS = regionprops(CC,'Area','BoundingBox','Centroid','PixelList','PixelIdxList');
% STATS = regionprops(CC,'Area','Centroid','PixelList');

%fprintf('Completed and found %d\n',CC.NumObjects);
% make into column vector

% preallocate other stats of interest
%%
if CC.NumObjects <= intmax('uint8')
    dataType = 'uint8';
elseif CC.NumObjects <= intmax('uint16')
    dataType = 'uint16';
elseif CC.NumObjects <= intmax('uint32')
    dataType = 'uint32';
else
    dataType = 'double';
end
% make new directory for output of surfaces
i_body = 1;
% I = I(:);
length_I = n_slices*skeyn*pkeyn;
clear I;

 % need to automate this selection
pkeys = job_meta.pkey_min:job_meta.pkey_inc:job_meta.pkey_max;
skeys = job_meta.skey_min:job_meta.skey_inc:job_meta.skey_max;
NumObjects = CC.NumObjects;
clear CC
% make directory for storing results
% str_date = date;
% str_date = regexprep(str_date, '-', '');
% out_dir_top = strcat(job_meta.output_dir,'/horizons_top',str_date,'_',num2str(anomalous_threshold),'_',num2str(area_thresh_lt),'_',num2str(area_thresh_gt));
% out_dir_bot = strcat(job_meta.output_dir,'/horizons_bot',str_date,'_',num2str(anomalous_threshold),'_',num2str(area_thresh_lt),'_',num2str(area_thresh_gt));
% out_dir_img = strcat(job_meta.output_dir,'/horizons_img',str_date,'_',num2str(anomalous_threshold),'_',num2str(area_thresh_lt),'_',num2str(area_thresh_gt));
% mkdir(out_dir_top);
% mkdir(out_dir_bot);
% mkdir(out_dir_img);
% for ii = 1:1:NumObjects
%     %I_vol(CC.PixelIdxList{ii}) = STATS(ii).Area; % volume
%     
%     if STATS(ii).Area > area_thresh_gt && STATS(ii).Area < area_thresh_lt
%         fprintf('Creating surfaces for body %d of %d\n',ii,CC.NumObjects);
%         z_col = 2;
%         skey_col = 1;
%         pkey_col = 3;
%         
%         pos = STATS(ii).PixelList;
%         sort_z_as = sortrows(pos,z_col);
%         sort_z_ds = sortrows(pos,-z_col);
%         pos_uniq = unique(pos(:,[pkey_col skey_col]),'rows');
%         % col 1 skey
%         % col 2 z
%         % col 3 pekyn
%         pos_body_surface_top(:,1) = pkeys(pos_uniq(:,1));
%         pos_body_surface_top(:,2) = skeys(pos_uniq(:,2));
%         
%         [~,Locb] = ismember(pos_uniq,sort_z_ds(:,[pkey_col skey_col]),'rows');
%         pos_body_surface_top(:,3) = (sort_z_ds(Locb,2)-1).*job_meta.s_rate/1000;
%         
%         pos_body_surface_bot(:,1) = pkeys(pos_uniq(:,1));
%         pos_body_surface_bot(:,2) = skeys(pos_uniq(:,2));
%         [~,Locb] = ismember(pos_uniq,sort_z_as(:,[3 1]),'rows');
%         pos_body_surface_bot(:,3) = (sort_z_as(Locb,2)-1).*job_meta.s_rate/1000;
%         % add mean depth to filename
%         top_asc = strcat(out_dir_top,'/anom_surface_body_',num2str(ii),'_top.txt');
%         bot_asc = strcat(out_dir_bot,'/anom_surface_body_',num2str(ii),'_bot.txt');
%         
%         area = repmat(STATS(ii).Area,size(pos_body_surface_bot(:,3),1),1); 
%         id = repmat(ii,size(pos_body_surface_bot(:,3),1),1); 
%         
%         dlmwrite(top_asc,[pos_body_surface_top area id],'delimiter','\t','precision', 8)
%         dlmwrite(bot_asc,[pos_body_surface_bot area id],'delimiter','\t','precision', 8)
%         
% %         h1a = scatter(pos_body_surface_top(:,1),pos_body_surface_top(:,2),6,pos_body_surface_top(:,3));
% %         %hold all
% %         %contour(x,y,z,'--w','LineWidth',2) 
% %         %hold off
% %         % add STATS(ii).Centroid(1)
% %         title(sprintf('Body: %d,\nAnomalous threshold: %f, \n Voxel cut-off: %f \n Centroid depth: %d\n' ...
% %             ,ii,anomalous_threshold,area_thresh,STATS(ii).Centroid(2)*job_meta.s_rate/1000),'interpreter','none');
% %         xlabel('Inline')
% %         ylabel('Crossline')
% %         colormap(flipud(jet))
% %         %c = colorbar;
% %         %ylabel(c,'Z') 
% %         %text(6500,6500,report_header)
% %         contour_fig = strcat(out_dir_img,'/body_',num2str(ii),'_',num2str(STATS(ii).Centroid(3)),'_',num2str(STATS(ii).Centroid(1)),'.eps'); 
% %         saveas(h1a,contour_fig,'eps') 
%         
%         % write out stats in filename and file list
% %         top_bulk = strcat(job_meta.output_dir,'/surface_body_bulk_',num2str(anomalous_threshold),'_top.txt');
% %         bot_bulk = strcat(job_meta.output_dir,'/surface_body_bulk_',num2str(anomalous_threshold),'_bot.txt');
% %         
% %         dlmwrite(top_bulk, [top_bulk; pos_body_surface_top],'delimiter','\t','precision', 8,'-append') 
% %         dlmwrite(bot_bulk, [bot_bulk; pos_body_surface_bot],'delimiter','\t','precision', 8,'-append') 
%         
%         i_body = i_body + 1;
%         
%         clear pos sort_z_as sort_z_ds pos_uniq Locb pos_body_surface_top pos_body_surface_bot
%     end
% end
% I = reshape(I,n_slices,skeyn,pkeyn);
% make a volume of zeros
n_lin_vols = 2;

if length_I/n_lin_vols == floor(length_I/n_lin_vols);
    length_I_vol(1) = length_I/n_lin_vols;
    length_I_vol(2) = length_I-length_I_vol(1);
else
    
end

for i_vol = 1:size(length_I_vol,2)
    I_vol = zeros(length_I_vol(i_vol),1,dataType);
    for ii = 1:1:NumObjects
        if max(STATS(ii).PixelIdxList) <= length_I_vol(i_vol)
            I_vol(STATS(ii).PixelIdxList) = STATS(ii).Area;
            STATS(ii) = [];
        else
            IdxList(i_lin).PixelIdxList = STATS(ii).PixelIdxList - length_I_vol(i_vol);
            IdxList(i_lin).Area = STATS(ii).Area;
            STATS(ii) = [];
            i_lin = i_lin + 1;
        end
    end
end

I_vol = zeros(length_I,1,dataType);
for ii = 1:1:NumObjects
    if max(STATS(ii).PixelIdxList) < 
    I_vol(STATS(ii).PixelIdxList) = STATS(ii).Area; % volume
end

start_il = job_meta.pkey_min;
start_xl = job_meta.skey_min;
n_xls = skeyn;
n_samps = n_slices;

fid = fopen(strcat(job_meta.output_dir,sprintf('%s_anomalies_result_%s_sil%s_sxl%s_nxl%s_samp%s.bin',num2str(anomalous_threshold),'I_vol',num2str(start_il),num2str(start_xl),num2str(n_xls),num2str(n_samps))),'w');
fwrite(fid,I_vol,'float32');
fclose(fid);
clear I_vol

%%
I_id = zeros(length_I,1,dataType);
for ii = 1:1:NumObjects
    % anomaly depth
    I_id(STATS(ii).PixelIdxList) = ii;
end
% I_id = labelmatrix(CC);

% Write out ID volume
fid = fopen(strcat(job_meta.output_dir,sprintf('%s_anomalies_result_%s_sil%s_sxl%s_nxl%s_samp%s.bin',num2str(anomalous_threshold),'I_id',num2str(start_il),num2str(start_xl),num2str(n_xls),num2str(n_samps))),'w');
fwrite(fid,I_id,'float32');
fclose(fid);
clear I_id

%%
% I_depth = zeros(length_I,1,dataType);
% for ii = 1:1:CC.NumObjects
%     % anomaly depth
%     I_depth(CC.PixelIdxList{ii}) = (STATS(ii).Centroid(2)+...
%         n_slices)*(job_meta.s_rate/1000);
% end
% fid = fopen(strcat(job_meta.output_dir,sprintf('%s_anomalies_result2_%s_sil%s_sxl%s_nxl%s_samp%s.bin',num2str(anomalous_threshold),'I_depth',num2str(start_il),num2str(start_xl),num2str(n_xls),num2str(n_samps))),'w');
% fwrite(fid,I_depth,'float32');
% fclose(fid);
% clear I_depth
% 
% % I_azimuth_axis = zeros(length_I,1,dataType);
% % I_azimuth_vals = zeros(length_I,1,dataType);
% % for ii = 1:1:CC.NumObjects
% %     % Calculate the azimuth (brace yourself!)
% %     % sort by x
% %     pos_sort_x = sortrows(unique(STATS(ii).PixelList,'rows'),[3 1 2]);
% %     pos_x_min = pos_sort_x(1,:);
% %     pos_x_max = pos_sort_x(end,:);
% %     % sort by y
% %     pos_sort_y = sortrows(unique(STATS(ii).PixelList,'rows'),[1 3 2]);
% %     pos_y_min = pos_sort_y(1,:);
% %     pos_y_max = pos_sort_y(end,:);
% %     
% %     diff_x = (pos_x_max(3)-pos_x_min(3))^2+(pos_x_max(1)-pos_x_min(1))^2;
% %     diff_y = (pos_y_max(3)-pos_y_min(3))^2+(pos_y_max(1)-pos_y_min(1))^2;
% %     
% %     if diff_x > diff_y
% %         minpos(1,:) = pos_x_min;
% %         maxpos(1,:) = pos_x_max;
% %     else
% %         minpos(1,:) = pos_y_min;
% %         maxpos(1,:) = pos_y_max;
% %     end
% %     % Calculate azimuth
% %     opp = maxpos(1) - minpos(1); % y dimension
% %     adj = maxpos(3) - minpos(3); % x dimension
% %     azimuth = atan(opp/adj);
% %     if isnan(azimuth)
% %         azimuth = 0;
% %         opp = 0;
% %         adj = 0;
% %     end
% %     
% %     I_azimuth_axis(CC.PixelIdxList{ii}) = abs(sin(azimuth)*opp);
% %     azimuth = 180/pi()*azimuth;
% %     if azimuth < 0;
% %         azimuth = abs(azimuth)+90;
% %     end
% %     I_azimuth_vals(CC.PixelIdxList{ii}) = azimuth;
% % end
% % fid1 = fopen(strcat(job_meta.output_dir,sprintf('%s_anomalies_result_%s.bin',num2str(anomalous_threshold),'I_azimuth_axis')),'w');
% % fwrite(fid1,I_azimuth_axis,'float32');
% % fclose(fid1);
% % clear I_azimuth_axis
% % 
% % fid2 = fopen(strcat(job_meta.output_dir,sprintf('%s_anomalies_result_%s.bin',num2str(anomalous_threshold),'I_azimuth_vals')),'w');
% % fwrite(fid2,I_azimuth_vals,'float32');
% % fclose(fid2);
% % clear I_azimuth_vals
% 
% %
% % % if only 1 block save result
% %     % if process_files.n_blocks == 1
% %         results_mat = strcat(job_meta.output_dir,...
% %                 'connectivity_results_block_.mat');
% %         struct_out.results_out{1,1} = 'body_id';
% %         struct_out.results_out{2,1} = 'body_volume';
% %         struct_out.results_out{3,1} = 'body_depth';
% %         struct_out.results_out{4,1} = 'body_azimuth';
% %         struct_out.results_out{1,2} = I_id;
% %         struct_out.results_out{2,2} = I_vol;
% %         struct_out.results_out{3,2} = I_depth;
% %         struct_out.results_out{4,2} = I_azimuth_vals;
% %         struct_out.anomalous_threshold = anomalous_threshold;
% %         save(results_mat,'-struct','struct_out','-v7.3');
% %
% %         for i_res = 1:1:4
% %             %output_dir = '/data/TZA/dtect/TZA_Kusini_Outboard_10x10/matlab_out/single_block/';
% %             fid = fopen(strcat(job_meta.output_dir,sprintf('%s_anomalies_result_%s.bin',num2str(anomalous_threshold),struct_out.results_out{i_res,1})),'w');
% %             fwrite(fid,struct_out.results_out{i_res,2},'float32');
% %             fclose(fid);
% %         end
% 
% % should just save as binary or segy
% %     else % multiple blocks
% %         I_vol = reshape(I_vol,n_slices,nil*nxl);
% %         I_id = reshape(I_id,n_slices,nil*nxl);
% %         I_azimuth_axis = reshape(I_azimuth_axis,n_slices,nil*nxl);
% %         I_azimuth_vals = reshape(I_azimuth_vals,n_slices,nil*nxl);
% %         I_depth = reshape(I_depth,n_slices,nil*nxl);
% %
% %         % save results from each block as matlab .mat binaries
% %         temp_results_mat = strcat(process_files.path_for_blocks,'.temp_',...
% %             process_files.func_name,'_block_',num2str(process_positions.block_id),'.mat');
% %         save(temp_results_mat,'anomalous_threshold','I_vol','I_id','I_depth','I_azimuth_axis','I_azimuth_vals','-v7.3');
% %
% %         temp_results_mat2 = strcat(process_files.path_for_blocks,'temp_',...
% %             process_files.func_name,'_block_',num2str(process_positions.block_id),'.mat');
% %         system_for = sprintf('mv %s %s',temp_results_mat,temp_results_mat2);
% %         system(system_for);
% %
% %         if join_block_id == process_positions.block_id; % this will be the node that joins things together
% %             clearvars -except process_files
% %             merge_anomalous_body_connector(process_files.func_name,process_files.path_for_blocks);
% %         end
% %     end
% 
% % ls *eps > eps_files.txt
% % gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -dGraphicsAlphaBits=4 -sOutputFile=body_30426482_5167.4298_7005.2478.png body_30426482_5167.4298_7005.2478.eps
% 
% % awk '{print "gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -dGraphicsAlphaBits=4 -sOutputFile="$0".png " $0}' eps_files.txt > convert_files.txt
% % Then run convert_files.txt.







% STATS = regionprops(CC,'Area','Centroid','PixelList');

%fprintf('Completed and found %d\n',CC.NumObjects);
% make into column vector

% preallocate other stats of interest
%%
if CC.NumObjects <= intmax('uint8')
    dataType = 'uint8';
elseif CC.NumObjects <= intmax('uint16')
    dataType = 'uint16';
elseif CC.NumObjects <= intmax('uint32')
    dataType = 'uint32';
else
    dataType = 'double';
end
% make new directory for output of surfaces
i_body = 1;
% I = I(:);
length_I = n_slices*skeyn*pkeyn;
clear I;

 % need to automate this selection
pkeys = job_meta.pkey_min:job_meta.pkey_inc:job_meta.pkey_max;
skeys = job_meta.skey_min:job_meta.skey_inc:job_meta.skey_max;
NumObjects = CC.NumObjects;
clear CC
% make directory for storing results
% str_date = date;
% str_date = regexprep(str_date, '-', '');
% out_dir_top = strcat(job_meta.output_dir,'/horizons_top',str_date,'_',num2str(anomalous_threshold),'_',num2str(area_thresh_lt),'_',num2str(area_thresh_gt));
% out_dir_bot = strcat(job_meta.output_dir,'/horizons_bot',str_date,'_',num2str(anomalous_threshold),'_',num2str(area_thresh_lt),'_',num2str(area_thresh_gt));
% out_dir_img = strcat(job_meta.output_dir,'/horizons_img',str_date,'_',num2str(anomalous_threshold),'_',num2str(area_thresh_lt),'_',num2str(area_thresh_gt));
% mkdir(out_dir_top);
% mkdir(out_dir_bot);
% mkdir(out_dir_img);
% for ii = 1:1:NumObjects
%     %I_vol(CC.PixelIdxList{ii}) = STATS(ii).Area; % volume
%     
%     if STATS(ii).Area > area_thresh_gt && STATS(ii).Area < area_thresh_lt
%         fprintf('Creating surfaces for body %d of %d\n',ii,CC.NumObjects);
%         z_col = 2;
%         skey_col = 1;
%         pkey_col = 3;
%         
%         pos = STATS(ii).PixelList;
%         sort_z_as = sortrows(pos,z_col);
%         sort_z_ds = sortrows(pos,-z_col);
%         pos_uniq = unique(pos(:,[pkey_col skey_col]),'rows');
%         % col 1 skey
%         % col 2 z
%         % col 3 pekyn
%         pos_body_surface_top(:,1) = pkeys(pos_uniq(:,1));
%         pos_body_surface_top(:,2) = skeys(pos_uniq(:,2));
%         
%         [~,Locb] = ismember(pos_uniq,sort_z_ds(:,[pkey_col skey_col]),'rows');
%         pos_body_surface_top(:,3) = (sort_z_ds(Locb,2)-1).*job_meta.s_rate/1000;
%         
%         pos_body_surface_bot(:,1) = pkeys(pos_uniq(:,1));
%         pos_body_surface_bot(:,2) = skeys(pos_uniq(:,2));
%         [~,Locb] = ismember(pos_uniq,sort_z_as(:,[3 1]),'rows');
%         pos_body_surface_bot(:,3) = (sort_z_as(Locb,2)-1).*job_meta.s_rate/1000;
%         % add mean depth to filename
%         top_asc = strcat(out_dir_top,'/anom_surface_body_',num2str(ii),'_top.txt');
%         bot_asc = strcat(out_dir_bot,'/anom_surface_body_',num2str(ii),'_bot.txt');
%         
%         area = repmat(STATS(ii).Area,size(pos_body_surface_bot(:,3),1),1); 
%         id = repmat(ii,size(pos_body_surface_bot(:,3),1),1); 
%         
%         dlmwrite(top_asc,[pos_body_surface_top area id],'delimiter','\t','precision', 8)
%         dlmwrite(bot_asc,[pos_body_surface_bot area id],'delimiter','\t','precision', 8)
%         
% %         h1a = scatter(pos_body_surface_top(:,1),pos_body_surface_top(:,2),6,pos_body_surface_top(:,3));
% %         %hold all
% %         %contour(x,y,z,'--w','LineWidth',2) 
% %         %hold off
% %         % add STATS(ii).Centroid(1)
% %         title(sprintf('Body: %d,\nAnomalous threshold: %f, \n Voxel cut-off: %f \n Centroid depth: %d\n' ...
% %             ,ii,anomalous_threshold,area_thresh,STATS(ii).Centroid(2)*job_meta.s_rate/1000),'interpreter','none');
% %         xlabel('Inline')
% %         ylabel('Crossline')
% %         colormap(flipud(jet))
% %         %c = colorbar;
% %         %ylabel(c,'Z') 
% %         %text(6500,6500,report_header)
% %         contour_fig = strcat(out_dir_img,'/body_',num2str(ii),'_',num2str(STATS(ii).Centroid(3)),'_',num2str(STATS(ii).Centroid(1)),'.eps'); 
% %         saveas(h1a,contour_fig,'eps') 
%         
%         % write out stats in filename and file list
% %         top_bulk = strcat(job_meta.output_dir,'/surface_body_bulk_',num2str(anomalous_threshold),'_top.txt');
% %         bot_bulk = strcat(job_meta.output_dir,'/surface_body_bulk_',num2str(anomalous_threshold),'_bot.txt');
% %         
% %         dlmwrite(top_bulk, [top_bulk; pos_body_surface_top],'delimiter','\t','precision', 8,'-append') 
% %         dlmwrite(bot_bulk, [bot_bulk; pos_body_surface_bot],'delimiter','\t','precision', 8,'-append') 
%         
%         i_body = i_body + 1;
%         
%         clear pos sort_z_as sort_z_ds pos_uniq Locb pos_body_surface_top pos_body_surface_bot
%     end
% end
% I = reshape(I,n_slices,skeyn,pkeyn);
% make a volume of zeros
n_lin_vols = 2;

if length_I/n_lin_vols == floor(length_I/n_lin_vols);
    length_I_vol(1) = length_I/n_lin_vols;
    length_I_vol(2) = length_I-length_I_vol(1);
else
    
end

for i_vol = 1:size(length_I_vol,2)
    I_vol = zeros(length_I_vol(i_vol),1,dataType);
    for ii = 1:1:NumObjects
        if max(STATS(ii).PixelIdxList) <= length_I_vol(i_vol)
            I_vol(STATS(ii).PixelIdxList) = STATS(ii).Area;
            STATS(ii) = [];
        else
            IdxList(i_lin).PixelIdxList = STATS(ii).PixelIdxList - length_I_vol(i_vol);
            IdxList(i_lin).Area = STATS(ii).Area;
            STATS(ii) = [];
            i_lin = i_lin + 1;
        end
    end
end

I_vol = zeros(length_I,1,dataType);
for ii = 1:1:NumObjects
    if max(STATS(ii).PixelIdxList) < 
    I_vol(STATS(ii).PixelIdxList) = STATS(ii).Area; % volume
end

start_il = job_meta.pkey_min;
start_xl = job_meta.skey_min;
n_xls = skeyn;
n_samps = n_slices;

fid = fopen(strcat(job_meta.output_dir,sprintf('%s_anomalies_result_%s_sil%s_sxl%s_nxl%s_samp%s.bin',num2str(anomalous_threshold),'I_vol',num2str(start_il),num2str(start_xl),num2str(n_xls),num2str(n_samps))),'w');
fwrite(fid,I_vol,'float32');
fclose(fid);
clear I_vol

%%
I_id = zeros(length_I,1,dataType);
for ii = 1:1:NumObjects
    % anomaly depth
    I_id(STATS(ii).PixelIdxList) = ii;
end
% I_id = labelmatrix(CC);

% Write out ID volume
fid = fopen(strcat(job_meta.output_dir,sprintf('%s_anomalies_result_%s_sil%s_sxl%s_nxl%s_samp%s.bin',num2str(anomalous_threshold),'I_id',num2str(start_il),num2str(start_xl),num2str(n_xls),num2str(n_samps))),'w');
fwrite(fid,I_id,'float32');
fclose(fid);
clear I_id

%%
% I_depth = zeros(length_I,1,dataType);
% for ii = 1:1:CC.NumObjects
%     % anomaly depth
%     I_depth(CC.PixelIdxList{ii}) = (STATS(ii).Centroid(2)+...
%         n_slices)*(job_meta.s_rate/1000);
% end
% fid = fopen(strcat(job_meta.output_dir,sprintf('%s_anomalies_result2_%s_sil%s_sxl%s_nxl%s_samp%s.bin',num2str(anomalous_threshold),'I_depth',num2str(start_il),num2str(start_xl),num2str(n_xls),num2str(n_samps))),'w');
% fwrite(fid,I_depth,'float32');
% fclose(fid);
% clear I_depth
% 
% % I_azimuth_axis = zeros(length_I,1,dataType);
% % I_azimuth_vals = zeros(length_I,1,dataType);
% % for ii = 1:1:CC.NumObjects
% %     % Calculate the azimuth (brace yourself!)
% %     % sort by x
% %     pos_sort_x = sortrows(unique(STATS(ii).PixelList,'rows'),[3 1 2]);
% %     pos_x_min = pos_sort_x(1,:);
% %     pos_x_max = pos_sort_x(end,:);
% %     % sort by y
% %     pos_sort_y = sortrows(unique(STATS(ii).PixelList,'rows'),[1 3 2]);
% %     pos_y_min = pos_sort_y(1,:);
% %     pos_y_max = pos_sort_y(end,:);
% %     
% %     diff_x = (pos_x_max(3)-pos_x_min(3))^2+(pos_x_max(1)-pos_x_min(1))^2;
% %     diff_y = (pos_y_max(3)-pos_y_min(3))^2+(pos_y_max(1)-pos_y_min(1))^2;
% %     
% %     if diff_x > diff_y
% %         minpos(1,:) = pos_x_min;
% %         maxpos(1,:) = pos_x_max;
% %     else
% %         minpos(1,:) = pos_y_min;
% %         maxpos(1,:) = pos_y_max;
% %     end
% %     % Calculate azimuth
% %     opp = maxpos(1) - minpos(1); % y dimension
% %     adj = maxpos(3) - minpos(3); % x dimension
% %     azimuth = atan(opp/adj);
% %     if isnan(azimuth)
% %         azimuth = 0;
% %         opp = 0;
% %         adj = 0;
% %     end
% %     
% %     I_azimuth_axis(CC.PixelIdxList{ii}) = abs(sin(azimuth)*opp);
% %     azimuth = 180/pi()*azimuth;
% %     if azimuth < 0;
% %         azimuth = abs(azimuth)+90;
% %     end
% %     I_azimuth_vals(CC.PixelIdxList{ii}) = azimuth;
% % end
% % fid1 = fopen(strcat(job_meta.output_dir,sprintf('%s_anomalies_result_%s.bin',num2str(anomalous_threshold),'I_azimuth_axis')),'w');
% % fwrite(fid1,I_azimuth_axis,'float32');
% % fclose(fid1);
% % clear I_azimuth_axis
% % 
% % fid2 = fopen(strcat(job_meta.output_dir,sprintf('%s_anomalies_result_%s.bin',num2str(anomalous_threshold),'I_azimuth_vals')),'w');
% % fwrite(fid2,I_azimuth_vals,'float32');
% % fclose(fid2);
% % clear I_azimuth_vals
% 
% %
% % % if only 1 block save result
% %     % if process_files.n_blocks == 1
% %         results_mat = strcat(job_meta.output_dir,...
% %                 'connectivity_results_block_.mat');
% %         struct_out.results_out{1,1} = 'body_id';
% %         struct_out.results_out{2,1} = 'body_volume';
% %         struct_out.results_out{3,1} = 'body_depth';
% %         struct_out.results_out{4,1} = 'body_azimuth';
% %         struct_out.results_out{1,2} = I_id;
% %         struct_out.results_out{2,2} = I_vol;
% %         struct_out.results_out{3,2} = I_depth;
% %         struct_out.results_out{4,2} = I_azimuth_vals;
% %         struct_out.anomalous_threshold = anomalous_threshold;
% %         save(results_mat,'-struct','struct_out','-v7.3');
% %
% %         for i_res = 1:1:4
% %             %output_dir = '/data/TZA/dtect/TZA_Kusini_Outboard_10x10/matlab_out/single_block/';
% %             fid = fopen(strcat(job_meta.output_dir,sprintf('%s_anomalies_result_%s.bin',num2str(anomalous_threshold),struct_out.results_out{i_res,1})),'w');
% %             fwrite(fid,struct_out.results_out{i_res,2},'float32');
% %             fclose(fid);
% %         end
% 
% % should just save as binary or segy
% %     else % multiple blocks
% %         I_vol = reshape(I_vol,n_slices,nil*nxl);
% %         I_id = reshape(I_id,n_slices,nil*nxl);
% %         I_azimuth_axis = reshape(I_azimuth_axis,n_slices,nil*nxl);
% %         I_azimuth_vals = reshape(I_azimuth_vals,n_slices,nil*nxl);
% %         I_depth = reshape(I_depth,n_slices,nil*nxl);
% %
% %         % save results from each block as matlab .mat binaries
% %         temp_results_mat = strcat(process_files.path_for_blocks,'.temp_',...
% %             process_files.func_name,'_block_',num2str(process_positions.block_id),'.mat');
% %         save(temp_results_mat,'anomalous_threshold','I_vol','I_id','I_depth','I_azimuth_axis','I_azimuth_vals','-v7.3');
% %
% %         temp_results_mat2 = strcat(process_files.path_for_blocks,'temp_',...
% %             process_files.func_name,'_block_',num2str(process_positions.block_id),'.mat');
% %         system_for = sprintf('mv %s %s',temp_results_mat,temp_results_mat2);
% %         system(system_for);
% %
% %         if join_block_id == process_positions.block_id; % this will be the node that joins things together
% %             clearvars -except process_files
% %             merge_anomalous_body_connector(process_files.func_name,process_files.path_for_blocks);
% %         end
% %     end
% 
% % ls *eps > eps_files.txt
% % gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -dGraphicsAlphaBits=4 -sOutputFile=body_30426482_5167.4298_7005.2478.png body_30426482_5167.4298_7005.2478.eps
% 
% % awk '{print "gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -dGraphicsAlphaBits=4 -sOutputFile="$0".png " $0}' eps_files.txt > convert_files.txt
% % Then run convert_files.txt.


end
end
