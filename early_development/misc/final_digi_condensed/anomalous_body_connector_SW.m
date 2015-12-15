%function anomalous_body_connector_SW(block_mat_all,block_mat,process_files_mat,anomalous_threshold,connectivity,join_block_id)
function anomalous_body_connector_SW(seismic_mat_path,i_block,n_blocks,binary_slice_path,anomalous_threshold,connectivity,join_block_id)
%block_mat - information about the traces to process in this block
% process_files_mat - structures for the files to be read in
% anomalous_threshold - between 0 and 1
% connectivity - 6, 18 or 26
% join_block_id - the id of the block that will combine results, set to 0
% if no

anomalous_threshold = str2num(anomalous_threshold);
connectivity = str2num(connectivity);
join_block_id = str2num(join_block_id);
i_block = str2num(i_block);
n_blocks = str2num(n_blocks);


% make the processing positions for the function
% process_files = load(process_files_mat,'-mat');
%  seismic = load(seismic_mat_path);
tmp_seismic = fread(fopen(seismic_mat_path,'r'),'double');
seismic.n_samples = tmp_seismic(1,1);
seismic.trace_ilxl_bytes = reshape(tmp_seismic(2:end),3,[])';
seismic.n_traces = size(seismic.trace_ilxl_bytes,1);
nil = length(unique(seismic.trace_ilxl_bytes(:,1)));
nxl = length(unique(seismic.trace_ilxl_bytes(:,1)));

% Write info to log file
% log_fid = fopen(strcat(seismic.output_dir,'psalm_log.txt'),'a');
% fprintf(log_fid,sprintf('\nAnomalous Body Connector initiated: %s.\n',datestr(clock)));
% fprintf(log_fid,sprintf('Anomalous Threshold: %d.\n',anomalous_threshold));
% fprintf(log_fid,sprintf('Connectivity: %d. \n',connectivity));
% fprintf(log_fid,sprintf('Join_block_id: %d. \n\n\',join_block_id));
% fclose(log_fid);

% Read - which should be slice order from the anomaly calculator function
%window_length = 0;
% file_to_open = sprintf('%ssas_result_slice_%d.bin',binary_slice_path,i_block);
file_to_open = sprintf('/data/Global/dtect/SRW_test_dataset/SEGY/matlab/test_n_blocks_100/digi_results/minimum_energy_slices/digi_minimum_energy_eer_projection_slices_block_%d.bin',i_block);
fid = fopen(file_to_open,'r');
ns = fread(fid,1,'float32');
fclose(fid);

% Initialize volume to read data into
slice_vol = false(nil,nxl,ns);

% Section to read data
%[slices] = node_binary_read_slices_lite(seismic,i_block,n_blocks,binary_slice_path,window_length,ns);

for ij = 1:1:n_blocks
    sas_slice_fid = fopen(sprintf('%ssas_result_slice_%d.bin',binary_slice_path,ij));
    slice = reshape(fread(sas_slice_fid,'float'),nil,nxl,[]);        
    fclose(sas_slice_fid);
    slice_z = size(slice,3);
    slice = slice >= anomalous_threshold;
    %     if ij == 1
    %         top = 1;
    %     else
    %         top = (ij-1)*slice_z+1;
    %     end
    %        base = ij*slice_z;
    if ij == 1
        top = 1;
        base = slice_z;
    else
        top = base + 1;
        base = top + slice_z - 1;
    end    
    slice_vol(:,:,top:base) = slice(:,:,:);
end

CC = bwconncomp(slice_vol,connectivity);

stat{1} = 'Area';
stat{2} = 'BoundingBox';
stat{3} = 'Centroid';
stat{4} = 'PixelList';
stat{5} = 'PixelIdxList';
L = labelmatrix(CC);

% for ii=1:1:ns    
%     slice = reshape(L(:,:,ii),801,801); imagesc(slice); pause(0.1);
%     orient = regionprops(slice,'Orientation');
% end

for jk = 1:5
    
    results_out{jk,1} = stat(jk);
    results_out{jk,2} = reshape(regionprops(CC,stat(jk)),nil*nxl,ns);
end

    test = 1;
    % Run connectivity algorithm:
    
    % [anom_traces_logic process_positions] = node_segy_read_traces(block_mat,process_files_mat,1);
    % [anom_traces_logic process_positions] = node_segy_read_traces(block_mat_all,block_mat,process_files_mat, 1); % (block_mat_all,block_mat,process_files_mat,index_file_read)
    %
    % if process_positions.ilxl_aperture ~= 0
    %     return;
    % else
    %     n_pos = (2*process_positions.ilxl_aperture+1)^2;
    % end
    
    % Position information
    % nil = length(unique(process_positions.ilxl_grid(:,1)))*sqrt(n_pos);
    % nxl = length(unique(process_positions.ilxl_grid(:,2)))*sqrt(n_pos);
    % nz = length(process_positions.z_pos);

% nil = seismic.n_iline;
% nxl = seismic.n_xline;
% nz = 50;

% 
% % Make input data in a 3D arrray for connectivity tool
% I = false(nz, ...
%     nxl, ...
%     nil);
% 
% 
% slices = slices >= anomalous_threshold; % create a logical array based on user threshold   
% 
% 
% 
% for ii = 1:1:nz
%     I(ii,:,:) = reshape(anom_traces_logic(ii,:),nxl,nil);
% end



% 
% % run connectivity and stats algorithm
% fprintf('-- Block %d --\n',process_positions.block_id);
% CC = bwconncomp(I,connectivity);
% % L = labelmatrix(CC);
% STATS = regionprops(CC,'Area','BoundingBox','Centroid','PixelList','PixelIdxList');
% 
% % make into column vector
% I = I(:);
% % preallocate other stats of interest
% I_vol = zeros(length(I),1);
% I_id = zeros(length(I),1);
% I_depth = zeros(length(I),1);
% I_azimuth_axis = zeros(length(I),1);
% I_azimuth_vals = zeros(length(I),1);
% 
% fprintf('-- Finished calculating connectivity for %d bodies --\n',CC.NumObjects);
% fprintf('-- Calculating volume, depth and azimuth --\n');
% 
% for ii = 1:1:CC.NumObjects
%     I_vol(CC.PixelIdxList{ii}) = STATS(ii).Area; % volume
%     
%     % generate random numbers for ids in each block
%     a = 1;
%     b = 500;
%     r = a + (b-a).*rand(1);
%     I_id(CC.PixelIdxList{ii}) = floor(ii*r);
%     
%     % anomaly depth
%     I_depth(CC.PixelIdxList{ii}) = (STATS(ii).Centroid(2)+...
%         (process_positions.block_id-1)*nz)*(process_files.s_rate/1000);
%     
%     % Calculate the azimuth (brace yourself!)
%     % sort by x
%         pos_sort_x = sortrows(unique(STATS(ii).PixelList,'rows'),[3 1 2]);
%         pos_x_min = pos_sort_x(1,:);
%         pos_x_max = pos_sort_x(end,:);
%         % sort by y
%         pos_sort_y = sortrows(unique(STATS(ii).PixelList,'rows'),[1 3 2]);
%         pos_y_min = pos_sort_y(1,:);
%         pos_y_max = pos_sort_y(end,:);
% 
%         diff_x = (pos_x_max(3)-pos_x_min(3))^2+(pos_x_max(1)-pos_x_min(1))^2;
%         diff_y = (pos_y_max(3)-pos_y_min(3))^2+(pos_y_max(1)-pos_y_min(1))^2;
% 
%         if diff_x > diff_y
%             minpos{ii}(1,:) = pos_x_min;
%             maxpos{ii}(1,:) = pos_x_max;
%         else
%             minpos{ii}(1,:) = pos_y_min;
%             maxpos{ii}(1,:) = pos_y_max;
%         end
%         % Calculate azimuth
%         opp = maxpos{ii}(1) - minpos{ii}(1); % y dimension    
%         adj = maxpos{ii}(3) - minpos{ii}(3); % x dimension
%         azimuth = atan(opp/adj);
%         if isnan(azimuth)
%            azimuth = 0; 
%            opp = 0;
%            adj = 0;
%         end
%         
%         I_azimuth_axis(CC.PixelIdxList{ii}) = abs(sin(azimuth)*opp); 
%         azimuth = 180/pi()*azimuth;
%         if azimuth < 0;
%            azimuth = abs(azimuth)+90; 
%         end
%         I_azimuth_vals(CC.PixelIdxList{ii}) = azimuth;
%     
%     if (floor(ii/5000) == ii/5000)
%         fprintf('%d%% complete \n',round((ii/CC.NumObjects)*100));
%     elseif ii == CC.NumObjects
%         fprintf('Completed block %d \n',process_positions.block_id);
%     end
% end
% 
% % if only 1 block save result
%     if process_files.n_blocks == 1
%         results_mat = strcat(process_files.path_for_blocks,...
%                 process_files.func_name,'_results_block_',...
%                 num2str(process_files.n_blocks),'.mat');
%         struct_out.results_out{1,1} = 'body_id';
%         struct_out.results_out{2,1} = 'body_volume';
%         struct_out.results_out{3,1} = 'body_depth';
%         struct_out.results_out{4,1} = 'body_azimuth';
%         struct_out.results_out{1,2} = I_id;
%         struct_out.results_out{2,2} = I_vol;
%         struct_out.results_out{3,2} = I_depth;
%         struct_out.results_out{4,2} = I_azimuth_vals;
%         struct_out.anomalous_threshold = anomalous_threshold;
%         save(results_mat,'-struct','struct_out','-v7.3');
%         
%         for i_res = 1:1:4
%             output_dir = '/data/TZA/dtect/TZA_Kusini_Outboard_10x10/matlab_out/single_block/';
%             fid = fopen(strcat(output_dir,sprintf('%s_anomalies_result_%s.bin',num2str(anomalous_threshold),struct_out.results_out{i_res,1})),'w');
%             fwrite(fid,struct_out.results_out{i_res,2},'float32');
%             fclose(fid);
%         end
%         
%         % should just save as binary or segy
%     else % multiple blocks
%         I_vol = reshape(I_vol,nz,nil*nxl);
%         I_id = reshape(I_id,nz,nil*nxl);
%         I_azimuth_axis = reshape(I_azimuth_axis,nz,nil*nxl);
%         I_azimuth_vals = reshape(I_azimuth_vals,nz,nil*nxl);
%         I_depth = reshape(I_depth,nz,nil*nxl);   
% 
%         % save results from each block as matlab .mat binaries
%         temp_results_mat = strcat(process_files.path_for_blocks,'.temp_',...
%             process_files.func_name,'_block_',num2str(process_positions.block_id),'.mat');
%         save(temp_results_mat,'anomalous_threshold','I_vol','I_id','I_depth','I_azimuth_axis','I_azimuth_vals','-v7.3');
%         
%         temp_results_mat2 = strcat(process_files.path_for_blocks,'temp_',...
%             process_files.func_name,'_block_',num2str(process_positions.block_id),'.mat');
%         system_for = sprintf('mv %s %s',temp_results_mat,temp_results_mat2);
%         system(system_for);
%         
%         if join_block_id == process_positions.block_id; % this will be the node that joins things together
%             clearvars -except process_files
%             merge_anomalous_body_connector(process_files.func_name,process_files.path_for_blocks); 
%         end
%     end
end
