function water_bottom_flatten(block_mat_all,block_mat,process_files_mat,join_block_id)
% Automatic water bottom picker
% Input seismic traces
% Output index of water bottom
% Should make it check that input is seismic data

    process_files = load(process_files_mat,'-mat');
    
    % make the processing positions for the function
    [traces process_positions] = node_segy_read_traces(block_mat_all,block_mat,process_files_mat,1);

    % A gain is applied to suppress all events inversely with depth
    % The water bottom should then be the brightest event
    exponent = 50;
    inverse_gain = repmat((size(traces,1):-1:1)'.^exponent,1,size(traces,2));    
        
    % Returns the index of the maximum value down the trace
    [max_val,index_wb] = max(traces.*inverse_gain);
    
    index_wb(~isnan(max_val)) = index_wb(~isnan(max_val)) - 10;
    
    positions = process_positions.ilxl_pos(:,1:2);
   
    index_wb_mat = strcat(process_files.path_for_blocks,...
            process_files.func_name,'_block_',...
            num2str(process_positions.block_id),'.mat');
    save(index_wb_mat,'index_wb','positions','-v7.3');   
    
    % join results together and add to function process_filesmat file
    if isequal(str2double(join_block_id),process_positions.block_id); 
        % this will be the node that joins things together
        n_mat = 0;
        while n_mat <= process_files.n_blocks           
            system_for = sprintf('ls -B %s%s_block* | wc -l',process_files.path_for_blocks,process_files.func_name);
            [~,n_mat] = system(system_for);  
            n_mat = str2double(n_mat);
            if n_mat == process_files.n_blocks % we have the files we need
                for i_block = 1:1:n_mat
                    index_wb_mat = strcat(process_files.path_for_blocks,...
                            process_files.func_name,'_block_',...
                            num2str(i_block),'.mat');

                    load(index_wb_mat);
                    index_wb_load{i_block} = index_wb;
                    positions_wb_load{i_block} = positions';
                end
                index_wb = [cell2mat(positions_wb_load)' cell2mat(index_wb_load)'];  
                n_samples_wb = process_files.n_samples - min(index_wb(:,3));
                % Save to process file of function user wants to run
                process_files_mat = strcat(process_files.path_for_blocks, ...
                    process_files.func_to_run,'_process_files.mat');
                save(process_files_mat,'index_wb','n_samples_wb','-append','-v7.3');
                n_mat = process_files.n_blocks + 1; % break while loop     
            end         
        end
    end
end    
