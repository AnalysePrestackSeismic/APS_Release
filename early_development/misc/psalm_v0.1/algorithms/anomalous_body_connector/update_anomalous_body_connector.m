function update_anomalous_body_connector(func_name,path_for_blocks,i_block)
 
    fprintf('-- Updating statistics for block %d --\n',i_block); 
    temp_results_mat = strcat(path_for_blocks,...
        'temp_updated_',func_name,'_block_',...
        num2str(i_block),'.mat');  
    load(temp_results_mat); 
    connect_mat = strcat(path_for_blocks,...
                'connect_results.mat');
    connect_out = load(connect_mat);

    for i_id = 1:1:length(connect_out.vol(:,1));
       positions_to_update = (results_out{1,2} == connect_out.vol(i_id,1));                           

       vol_val = sum(connect_out.vol(i_id,2:end));
       depth_val = mean(connect_out.depth(i_id,2:end));
       azimuth_val = max(connect_out.azi(i_id,2:end));

       results_out{2,2}(positions_to_update) = vol_val;
       results_out{3,2}(positions_to_update) = depth_val;
       results_out{4,2}(positions_to_update) = azimuth_val;   


        if (floor(i_id/10) == i_id/10)
            fprintf('%d%% complete \n',round((i_id/length(ids(:,1))*100)));
        elseif i_id == length(ids(:,1))
            fprintf('Completed block %d.\n',i_block);
        end

    end

%   node_segy_write_traces(results_out,i_block,output_dir);

    results_mat = strcat(process_files.path_for_blocks,...
        process_files.func_name,'_results_block_',...
        num2str(i_block),'.mat');
    save(results_mat,'results_out','anomalous_threshold','-v7.3');

        
end