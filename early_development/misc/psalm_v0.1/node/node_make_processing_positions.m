function process_positions = node_make_processing_positions(block_mat_all,block_mat)
% for very large files this may need to run in a distributed manner    
    process_positions_all_tmp = load(block_mat_all,'-mat');
    process_positions_tmp = load(block_mat,'-mat');
    
    FN = {fieldnames(process_positions_all_tmp) fieldnames(process_positions_tmp)};
    VAL = {struct2cell(process_positions_all_tmp) struct2cell(process_positions_tmp)};
    
    FN = cat(1,FN{:});
    VAL = cat(1,VAL{:});
    process_positions = cell2struct(VAL, FN);
    
    % find bytes for the ilxl positions
    [Lia, Locb] = ...
        ismember(process_positions.ilxl_grid(:,1:2),process_positions.trace_ilxl_bytes(:,1:2), 'rows');
    % get trace bytes
    process_positions.ilxl_grid(:,3) = NaN(length(process_positions.ilxl_grid(:,1)),1);
    process_positions.ilxl_grid(Lia,3) = process_positions.trace_ilxl_bytes(nonzeros(Locb),3);

    if process_positions.ilxl_aperture > 0
        n_pos = (2*process_positions.ilxl_aperture+1)^2;  
        positions_right = 0:process_positions.ilxl_aperture_step:process_positions.ilxl_aperture;
        positions_left = fliplr(-process_positions.ilxl_aperture_step:-process_positions.ilxl_aperture_step:-process_positions.ilxl_aperture);    
        pos_diff_xl = repmat([positions_left positions_right],sqrt(n_pos),1);
        pos_diff_il = pos_diff_xl';
        pos_diff_il =  pos_diff_il(:).*process_positions.il_inc; % turn into vector
        pos_diff_xl =  pos_diff_xl(:).*process_positions.xl_inc;
        rep_il_pos = bsxfun(@plus,process_positions.ilxl_grid(:,1), ...
            pos_diff_il');            
        rep_xl_pos = bsxfun(@plus,process_positions.ilxl_grid(:,2), ...
            pos_diff_xl');
        rep_il_pos = rep_il_pos';
        rep_xl_pos = rep_xl_pos';
        rep_il_pos = rep_il_pos(:);
        rep_xl_pos = rep_xl_pos(:);
        temp_positions = [rep_il_pos rep_xl_pos];
        [Lia2, Locb2] = ...
            ismember(temp_positions(:,1:2),process_positions.trace_ilxl_bytes(:,1:2), 'rows');
        % get trace bytes
        temp_positions(:,3) = NaN(length(temp_positions(:,1)),1);
        temp_positions(Lia2,3) = process_positions.trace_ilxl_bytes(nonzeros(Locb2),3);
        process_positions.ilxl_pos = temp_positions;  
    else
        process_positions.ilxl_pos = process_positions.ilxl_grid;            
    end        

    z_diff = -process_positions.z_aperture:1:process_positions.z_aperture;
    % we add the aperture matrix to every element in the matrix
    temp_process_z_pos = bsxfun(@plus,process_positions.z_grid, ...
     z_diff);
    temp_process_z_pos = temp_process_z_pos';
    temp_process_z_pos = temp_process_z_pos(:); 
    % this defines acceptable sample numbers, i.e. positive
    z_samples = (1:1:process_positions.n_samples)'; 
    %look up and make NaN negative and too large sample numbers
    [~, Locb3] = ...
    ismember(temp_process_z_pos,z_samples, 'rows');
    process_positions.z_pos = nonzeros(Locb3);

    % start byte will include the z pos and aperture
    process_positions.start_byte = min(process_positions.ilxl_pos(:,3))+(min(process_positions.z_pos)-1).*4;
end        
