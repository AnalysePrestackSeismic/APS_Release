function merge_anomalous_body_connector(func_name,path_for_blocks)

start_point = pwd;
cd(path_for_blocks);  

process_files_mat = strcat(path_for_blocks,func_name,'_process_files.mat');
process_files = load(process_files_mat,'-mat');

fprintf('\n-- Joining results from %d blocks --\n',process_files.n_blocks); 
           
n_mat = 0;
while n_mat <= process_files.n_blocks
    system_for = sprintf('ls -B %stemp_%s_block* | wc -l',path_for_blocks,func_name);
    [~,n_mat] = system(system_for);  
    n_mat = str2double(n_mat);

    if n_mat == process_files.n_blocks % we have the files we need
        for i_block = 1:1:n_mat-1
            fprintf('-- Assessing connectivity at boundary between block %d and block %d --\n',... 
                i_block,i_block+1); 
                    % load result for block
                    if i_block == 1                
                        block_connect_mat = strcat(process_files.path_for_blocks,'temp_',...
                            process_files.func_name,'_block_',num2str(i_block),'.mat');
                        %block_connect_mat = strcat(process_files.path_for_blocks,...
                        %    process_files.func_name,'_results_block_',num2str(i_block),'.mat');
                        
                        block_load{1} = load(block_connect_mat);
                    else
                        block_load{1} = block_load{2};
                        block_load{1}.I_id = reshape(block_load{1}.I_id,nz2,ntraces);
                        block_load{1}.I_vol = reshape(block_load{1}.I_vol,nz2,ntraces);
                        block_load{1}.I_depth = reshape(block_load{1}.I_depth,nz2,ntraces);
                        block_load{1}.I_azimuth_axis = reshape(block_load{1}.I_azimuth_axis,nz2,ntraces);
                        block_load{1}.I_azimuth_vals = reshape(block_load{1}.I_azimuth_vals,nz2,ntraces);
                    end

                    block_connect_mat = strcat(process_files.path_for_blocks,'temp_',...
                            process_files.func_name,'_block_',num2str(i_block+1),'.mat');
                        
                    %block_connect_mat = strcat(process_files.path_for_blocks,...
                    %        process_files.func_name,'_results_block_',num2str(i_block+1),'.mat');   
                    block_load{2} = load(block_connect_mat);
                    
                    nz1 = size(block_load{1}.I_id,1);
                    nz2 = size(block_load{2}.I_id,1);
                    ntraces = size(block_load{1}.I_id,2)*size(block_load{1}.I_id,3);

                    % bottom from {1}
                    bot_bound = block_load{1}.I_id(end,:);
                    top_bound = block_load{2}.I_id(1,:);
                    
                    bot_bound = bot_bound(:);
                    top_bound = top_bound(:);
                    
                    % Create slices for other statistics
                    vol_bot_bound = block_load{1}.I_vol(end,:);
                    vol_top_bound = block_load{2}.I_vol(1,:);
                    vol_bot_bound = vol_bot_bound(:);
                    vol_top_bound = vol_top_bound(:);
                    
                    depth_bot_bound = block_load{1}.I_depth(end,:);
                    depth_top_bound = block_load{2}.I_depth(1,:);
                    
                    depth_bot_bound = depth_bot_bound(:);
                    depth_top_bound = depth_top_bound(:);
                    
                    azimuth_bot_bound = block_load{1}.I_azimuth_vals(end,:);
                    azimuth_top_bound = block_load{2}.I_azimuth_vals(1,:);
                    
                    azimuth_bot_bound = azimuth_bot_bound(:);
                    azimuth_top_bound = azimuth_top_bound(:);
                    
                    % Find join locations
                    logical_test = (top_bound.*bot_bound > 0);
                    ids(logical_test,1:8) = [bot_bound(logical_test),...
                                            top_bound(logical_test),...
                                            vol_bot_bound(logical_test),...
                                            vol_top_bound(logical_test),...
                                            depth_bot_bound(logical_test),...
                                            depth_top_bound(logical_test),...
                                            azimuth_bot_bound(logical_test),...
                                            azimuth_top_bound(logical_test)...
                                            ];
                    ids = unique(ids,'rows');
                    ids = ids(2:end,:);
                    
                    fprintf('- Updating IDs for %d bodies\n',length(ids(:,1))); 
                    if i_block == 1;
                        col = 2;
                    else
                        col = col + 1;
                    end
                    
                    block_load{1}.I_id = block_load{1}.I_id(:);
                    block_load{2}.I_id = block_load{2}.I_id(:); 
                    
                    %
                    % bsxfun(@eq,
                    
                    
                    for i_id = 1:1:length(ids(:,1))
                        % Find index of ids that connect across the
                        % boundary in bottom block
                        positions_to_find = (block_load{2}.I_id == ids(i_id,2));
                        % Set the ids of these bodies to be the same as in
                        % top block
                        block_load{2}.I_id(positions_to_find) = ids(i_id,1);
                    
                        connect_out.vol(i_id,1) = ids(i_id,1);
                        connect_out.vol(i_id,col) = sum(unique(ids(ids == ids(i_id),3)));
                        connect_out.vol(i_id,col+1) = sum(unique(ids(ids == ids(i_id),4)));
                        
                        connect_out.depth(i_id,1) = ids(i_id,1);
                        connect_out.depth(i_id,col) = mean(unique(ids(ids == ids(i_id),5)));
                        connect_out.depth(i_id,col+1) = mean(unique(ids(ids == ids(i_id),6)));
                        
                        connect_out.azi(i_id,1) = ids(i_id,1);
                        connect_out.azi(i_id,col) = mode(unique(ids(ids == ids(i_id),7)));
                        connect_out.azi(i_id,col+1) = mode(unique(ids(ids == ids(i_id),8))); 
                        
                        if (floor(i_id/10) == i_id/10)
                            fprintf('%d%% complete \n',round((i_id/length(ids(:,1))*100)));
                        elseif i_id == length(ids(:,1))
                            fprintf('Completed block %d. Saving Results... \n',i_block);
                        end
                    end
                    
                    fprintf('- Saving results for boundary between block %d and block %d\n\n',... 
                        i_block,i_block+1);  
                    
                    % save block_load{1}
                    temp_results_mat = strcat(process_files.path_for_blocks,...
                            'temp_updated_',process_files.func_name,'_block_',...
                            num2str(i_block),'.mat');
                    struct_out.results_out{1,1} = 'body_id';
                    struct_out.results_out{2,1} = 'body_volume';
                    struct_out.results_out{3,1} = 'body_depth';
                    struct_out.results_out{4,1} = 'body_azimuth';
                    struct_out.results_out{1,2} = block_load{1}.I_id;
                    struct_out.results_out{2,2} = block_load{1}.I_vol;
                    struct_out.results_out{3,2} = block_load{1}.I_depth;
                    struct_out.results_out{4,2} = block_load{1}.I_azimuth_vals;
                    struct_out.anomalous_threshold = block_load{1}.anomalous_threshold;
                    save(temp_results_mat,'-struct','struct_out','-v7.3');
                    
                    if i_block == n_mat-1 % save final block
                        temp_results_mat = strcat(process_files.path_for_blocks,...
                            'temp_updated_',process_files.func_name,'_block_',...
                            num2str(i_block+1),'.mat');                           
                        struct_out.results_out{1,1} = 'body_id';
                        struct_out.results_out{2,1} = 'body_volume';
                        struct_out.results_out{3,1} = 'body_depth';
                        struct_out.results_out{4,1} = 'body_azimuth';
                        struct_out.results_out{1,2} = block_load{2}.I_id;
                        struct_out.results_out{2,2} = block_load{2}.I_vol;
                        struct_out.results_out{3,2} = block_load{2}.I_depth;
                        struct_out.results_out{4,2} = block_load{2}.I_azimuth_vals;
                        struct_out.anomalous_threshold = block_load{2}.anomalous_threshold;
                        save(temp_results_mat,'-struct','struct_out','-v7.3');
                    end
                    
                    clearvars ids
                     
        end
        
        % Save unique rows results
        % connect_out.vol = unique(connect_out.vol,'rows');
        % connect_out.depth = unique(connect_out.depth,'rows');
        % connect_out.azi = unique(connect_out.azi,'rows');     
      
        connect_mat = strcat(process_files.path_for_blocks,...
            'connect_results.mat');
        save(connect_mat,'-struct','connect_out','process_files','n_mat','start_point','-v7.3');    
        
    n_mat = process_files.n_blocks + 1; % break while loop     
    end
end
cd(start_point);
end