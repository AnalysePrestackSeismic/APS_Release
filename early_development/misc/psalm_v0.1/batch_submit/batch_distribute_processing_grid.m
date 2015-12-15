% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% PSALM - Parallel Seismic Analysis Leveraging MATLAB
% Date: 03 October 2012
% Authors: James Selvage and Jonathan Edgar
function batch_distribute_processing_grid(seismic,batch_structure)
   
    distribute_type = algorithm_to_run(batch_structure.func_name);
    
    if strcmp('either',distribute_type)
        distribute_type = input('Enter distribution type (''trace'' or ''slice''): ');
    end

    switch distribute_type 
        case 'trace'
            % Distribute for inline access
            n_pos = length(batch_structure.processing_grid.ilxl_grid);
        case 'slice'
            if batch_structure.in_flatten_water_bottom == 1
                shift_n_samples = batch_structure.n_samples-min(batch_structure.index_wb(:,3));
                [~,nan_index] = min(abs(batch_structure.processing_grid.z_grid - shift_n_samples));  
                batch_structure.processing_grid.z_grid = batch_structure.processing_grid.z_grid(1:nan_index);
            end
                n_pos = length(batch_structure.processing_grid.z_grid);
    end
    
    n_pos_per_block = floor(n_pos/batch_structure.n_blocks);
    n_pos_left_over = n_pos-(n_pos_per_block*batch_structure.n_blocks);
    % block index number
    block_i = 1:1:batch_structure.n_blocks;
    n_pos_on_block = repmat(n_pos_per_block,1,batch_structure.n_blocks);
    % add back remaining traces to last block
    n_pos_on_block(end) = n_pos_on_block(end)+n_pos_left_over;
    % calculate end index for each block
    end_index = block_i.*n_pos_per_block;         
    end_index(end) = end_index(end)+n_pos_left_over;
    start_index = end_index-(n_pos_on_block)+1;
    index_block = [block_i' start_index' end_index'];
    
%     n_pos_per_block = floor(n_pos/(batch_structure.n_blocks-1));
%     n_pos_left_over = n_pos-(n_pos_per_block*(batch_structure.n_blocks-1));
%     % block index number
%     block_i = 1:1:batch_structure.n_blocks;
%     n_pos_on_block = repmat(n_pos_per_block,1,batch_structure.n_blocks);
%     % add back remaining traces to last block
%     n_pos_on_block(end) = n_pos_left_over;
%     % calculate end index for each block
%     end_index = block_i.*n_pos_per_block;         
%     end_index(end) = n_pos;
%     start_index = end_index-(n_pos_on_block)+1;
%     index_block = [block_i' start_index' end_index']; 
    
    process_positions.z_step = batch_structure.processing_grid.z_step;
    process_positions.ilxl_step = batch_structure.processing_grid.ilxl_step;
    process_positions.ilxl_aperture = batch_structure.processing_grid.ilxl_aperture;
    process_positions.ilxl_aperture_step = batch_structure.processing_grid.ilxl_aperture_step;
    process_positions.z_aperture = batch_structure.processing_grid.z_aperture;

    process_positions.n_samples = seismic.n_samples;
    process_positions.s_rate = seismic.s_rate;
    process_positions.n_traces = seismic.n_traces;
    process_positions.trace_ilxl_bytes = seismic.trace_ilxl_bytes;
    process_positions.extra_bytes_data = seismic.extra_bytes_data;
    process_positions.min_iline = seismic.min_iline;
    process_positions.max_iline = seismic.max_iline;
    process_positions.min_xline = seismic.min_xline;
    process_positions.max_xline = seismic.max_xline;
    process_positions.n_iline = seismic.n_iline;
    process_positions.n_xline = seismic.n_xline;
    process_positions.il_inc = seismic.il_inc;
    process_positions.xl_inc = seismic.xl_inc;
    process_positions.distribute_type = distribute_type;
    process_positions.trc_head_length = seismic.trc_head_length;
    
    block_mat_all = strcat(batch_structure.path_for_blocks,batch_structure.func_name,'_positions_block_all.mat');
    save(block_mat_all,'-struct','process_positions','-v7.3');
    
    clearvars process_positions
       
    for ii = 1:1:batch_structure.n_blocks       
        switch distribute_type
            case 'trace' % for ilxl
                process_positions.ilxl_grid = batch_structure.processing_grid.ilxl_grid(index_block(ii,2):index_block(ii,3),:);
                process_positions.z_grid = batch_structure.processing_grid.z_grid;
            case 'slice' % for slices
                process_positions.ilxl_grid = batch_structure.processing_grid.ilxl_grid;
                if batch_structure.in_flatten_water_bottom == 1
                    process_positions.z_grid = batch_structure.processing_grid.z_grid(index_block(ii,2):index_block(ii,3),:)+min(batch_structure.index_wb(:,3));
                else
                    process_positions.z_grid = batch_structure.processing_grid.z_grid(index_block(ii,2):index_block(ii,3),:);
                end
        end
        
        process_positions.block_id = ii;
        
        block_mat = strcat(batch_structure.path_for_blocks,batch_structure.func_name,'_positions_block_',num2str(ii),'.mat');
 
        save(block_mat,'-struct','process_positions','-v7.3');

    end
    
end