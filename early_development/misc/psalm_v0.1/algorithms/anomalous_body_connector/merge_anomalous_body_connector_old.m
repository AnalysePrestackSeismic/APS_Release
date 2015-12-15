function merge_anomalous_body_connector(func_name,path_for_blocks)
%% Inputs:
% func_name = name of function that created the blocks to be joined e.g.
% anomalous_body_connector. Used to look up the _process_files.mat

% path_for_blocks = directory containing blocks to be joined.

%% Function:
start_point = pwd;
cd(path_for_blocks);  

process_files_mat = strcat(path_for_blocks,func_name,'_process_files.mat');
process_files = load(process_files_mat,'-mat');

fprintf('\n-- Joining results from %d blocks --\n',process_files.n_blocks); 
           
n_mat = 0;
while n_mat <= process_files.n_blocks
    system_for = sprintf('ls -B temp_%s_block* | wc -l',func_name);
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
                    block_load{2} = load(block_connect_mat);
                    
                    nz1 = size(block_load{1}.I_id,1); % number of samples in first block loaded
                    nz2 = size(block_load{2}.I_id,1); % number of samples in second block loaded
                    ntraces = size(block_load{1}.I_id,2)*size(block_load{1}.I_id,3);

                    % bottom from {1}
                    bot_bound = block_load{1}.I_id(end,:);
                    top_bound = block_load{2}.I_id(1,:);
                    
                    bot_bound = bot_bound(:);
                    top_bound = top_bound(:);
                    
                    logical_test = (top_bound.*bot_bound > 0);
                    ids(logical_test,1:2) = [bot_bound(logical_test) top_bound(logical_test)]; % Sets up logical test - gives 2xn matrix of zeros where n = number of ids that fit the logical test
                    ids = unique(ids,'rows'); % Returns unique ids
                    ids = ids(2:end,1:2); % First row is always zero from the previous step - not a valid id so is removed
                    % in ids column 1 is top block and column 2 is bottom
                    % block
                    
                    fprintf('- Updating IDs for %d bodies\n',length(ids(:,1))); 
                    if i_block == 1;
                        col = 2;
                    else
                        col = col + 1;
                    end
                    
                    % make all inputs into column vectors
                    block_load{1}.I_id = block_load{1}.I_id(:);
                    block_load{1}.I_vol = block_load{1}.I_vol(:);
                    block_load{1}.I_depth = block_load{1}.I_depth(:);
                    block_load{1}.I_azimuth_axis = block_load{1}.I_azimuth_axis(:);
                    block_load{1}.I_azimuth_vals = block_load{1}.I_azimuth_vals(:);

                    block_load{2}.I_id = block_load{2}.I_id(:);
                    block_load{2}.I_vol = block_load{2}.I_vol(:);
                    block_load{2}.I_depth = block_load{2}.I_depth(:);
                    block_load{2}.I_azimuth_axis = block_load{2}.I_azimuth_axis(:); 
                    block_load{2}.I_azimuth_vals = block_load{2}.I_azimuth_vals(:); 
                    
                    for i_id = 1:1:length(ids(:,1)) % can this be vectorised
                        positions_to_find = (block_load{2}.I_id == ids(i_id,2)); % find index locations for "wanted" ids in bottom block
                        block_load{2}.I_id(positions_to_find) = ids(i_id,1); % set ids of positions_to_find to ids from top block
   
                        % Check that ids were updated correctly
                        %block_1_check = reshape(block_load{1}.I_id,nz1,ntraces);
                        %block_2_check = reshape(block_load{2}.I_id,nz1,ntraces);
                        
                        connect_out.vol(i_id,1) = ids(i_id,1);
                        connect_out.depth(i_id,1) = ids(i_id,1);
                        connect_out.azimuth_axis(i_id,1) = ids(i_id,1);
                        connect_out.azimuth_vals(i_id,1) = ids(i_id,1);
                         
                        positions_to_update = (block_load{1}.I_id == ids(i_id,1));
                        connect_out.vol(i_id,col) = sum(unique(block_load{1}.I_vol(positions_to_update)));
                        %connect_out.depth(i_id,col) = mean(block_load{1}.I_depth(positions_to_update));
                        connect_out.depth(i_id,col) = mean(block_load{1}.I_depth(positions_to_update))+((nz1*(i_block-1))*process_files.s_rate/1000);
                        
                        
                        nan_azi = ~isnan(block_load{1}.I_azimuth_axis(positions_to_update));
                        connect_out.azimuth_axis(i_id,col) = mode(block_load{1}.I_azimuth_axis(nan_azi));                     
                        positions_to_update = (block_load{1}.I_azimuth_axis == connect_out.azimuth_axis(i_id,col));
                        connect_out.azimuth_vals(i_id,col) = mode(block_load{1}.I_azimuth_vals(positions_to_update));
                        
                        positions_to_update = (block_load{2}.I_id == ids(i_id,1));
                        connect_out.vol(i_id,col+1) = sum(unique(block_load{2}.I_vol(positions_to_update)));
                        %connect_out.depth(i_id,col+1) = mean(block_load{2}.I_depth(positions_to_update));
                        connect_out.depth(i_id,col+1) = mean(block_load{2}.I_depth(positions_to_update))+((nz2*(i_block))*process_files.s_rate/1000);
                        
                        
                        
                        nan_azi = ~isnan(block_load{2}.I_azimuth_axis(positions_to_update));
                        connect_out.azimuth_axis(i_id,col+1) = mode(block_load{2}.I_azimuth_axis(nan_azi));
                        positions_to_update = (block_load{2}.I_azimuth_axis == connect_out.azimuth_axis(i_id,col+1));
                        connect_out.azimuth_vals(i_id,col+1) = mode(block_load{2}.I_azimuth_vals(positions_to_update));
                        
                        if (floor(i_id/10) == i_id/10)
                            fprintf('%d%% complete \n',round((i_id/length(ids(:,1))*100)));
                        elseif i_id == length(ids(:,1))
                            fprintf('Completed block %d \n',i_block);
                        end
                    end
                    
                    % save block_load{1}
                    temp_results_mat = strcat(process_files.path_for_blocks,...
                            'temp_updated_',process_files.func_name,'_block_',...
                            num2str(i_block),'.mat');
                    struct_out.results_out{1,1} = 'body_id';
                    struct_out.results_out{2,1} = 'body_volume';
                    struct_out.results_out{3,1} = 'body_depth';
                    struct_out.results_out{4,1} = 'body_azimuth';
                    struct_out.results_out{1,2} = reshape(block_load{1}.I_id,nz1,ntraces);
                    struct_out.results_out{2,2} = reshape(block_load{1}.I_vol,nz1,ntraces);
                    struct_out.results_out{3,2} = reshape(block_load{1}.I_depth,nz1,ntraces);
                    struct_out.results_out{4,2} = reshape(block_load{1}.I_azimuth_vals,nz1,ntraces);
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
                        struct_out.results_out{1,2} = reshape(block_load{2}.I_id,nz2,ntraces);
                        struct_out.results_out{2,2} = reshape(block_load{2}.I_vol,nz2,ntraces);
                        struct_out.results_out{3,2} = reshape(block_load{2}.I_depth,nz2,ntraces);
                        struct_out.results_out{4,2} = reshape(block_load{2}.I_azimuth_vals,nz2,ntraces);
                        struct_out.anomalous_threshold = block_load{2}.anomalous_threshold;
                        save(temp_results_mat,'-struct','struct_out','-v7.3');
                    end
                    
                    % save joining information
                    connect_mat = strcat(process_files.path_for_blocks,...
                        'connect_results.mat');
                    save(connect_mat,'-struct','connect_out','-v7.3');

                    fprintf('- Saving results for boundary between block %d and block %d\n\n',... 
                        i_block,i_block+1);          
        end       
        
        clearvars -except process_files n_mat start_point

        % Update other stats
        for i_block = 1:1:n_mat    
            fprintf('-- Updating statistics for block %d --\n',i_block); 
            temp_results_mat = strcat(process_files.path_for_blocks,...
                'temp_updated_',process_files.func_name,'_block_',...
                num2str(i_block),'.mat');  
            load(temp_results_mat); 
            connect_mat = strcat(process_files.path_for_blocks,...
                        'connect_results.mat');
            connect_out = load(connect_mat);
       
            for i_id = 1:1:length(connect_out.vol(:,1));
               positions_to_update = (results_out{1,2} == connect_out.vol(i_id,1));                           

               vol_val = sum(connect_out.vol(i_id,2:end));
               depth_val = mean(connect_out.depth(i_id,2:end));
               [~,azimuth_id] = max(connect_out.azimuth_axis(i_id,2:end));
               azimuth_val = connect_out.azimuth_vals(i_id,azimuth_id+1);
               
               results_out{2,2}(positions_to_update) = vol_val;
               results_out{3,2}(positions_to_update) = depth_val;
               results_out{4,2}(positions_to_update) = azimuth_val;                           
            end 

            results_mat = strcat(process_files.path_for_blocks,...
                process_files.func_name,'_results_block_',...
                num2str(i_block),'.mat');
            save(results_mat,'results_out','anomalous_threshold','-v7.3');
        end
    n_mat = process_files.n_blocks + 1; % break while loop     
    end
end   
cd(start_point);
end