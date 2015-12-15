function [] = sas_to_segy_sw(seismic_mat_path,i_block,n_blocks,block_path,wb_path,output_dir)

i_block = str2double(i_block);
n_blocks = str2double(n_blocks);

tmp_seismic = fread(fopen(seismic_mat_path,'r'),'double');
seismic.n_samples = tmp_seismic(1,1);
%seismic.trace_ilxl_bytes = reshape(typecast(tmp_seismic(2:end),'int32'),3,[]');
seismic.trace_ilxl_bytes = reshape(tmp_seismic(2:end),3,[])';
seismic.n_traces = size(seismic.trace_ilxl_bytes,1);



for block = 1:1:n_blocks
    
    for slice = 1:1:n_blocks
        
        % Read in data from file
        file_to_open = sprintf('%ssas_result_slice_%d_trc_block_%d.bin',block_path,slice,block);
        
        fid = fopen(file_to_open,'r');
        
        % Calculate the start and end traces of this block
        
        if block == n_blocks
            end_trace_idx = seismic.n_traces;
            start_trace_idx = 1+((block-1)*floor(seismic.n_traces/n_blocks));
            
        else
            start_trace_idx = 1+((block-1)*floor(seismic.n_traces/n_blocks));
            end_trace_idx = block*floor(seismic.n_traces/n_blocks);
        end
        
        
        % end_trace_idx = slice*floor(seismic.n_traces/n_blocks);
        
        
        
        % Read datazz
        n_traces_to_read = end_trace_idx-start_trace_idx+1;
        traces{slice,1} = fread(fid,'float32');
        %slices{slice,block} = fread(fid,'float32');
        fclose(fid);
        
        
        
        % Arrange data into correct results out format to be written to
        % SEGY
        
        
        traces{slice,1} = reshape(traces{slice,1}(:,1),size(traces{slice,1}(:,1),1)/n_traces_to_read,n_traces_to_read);
    end
    sas_result_trc = vertcat(traces{:});
    
    
    %% Section to unflatten data
    
    water_bottom = sprintf('%swaterbottom_idx_block_%d.bin',wb_path,block);
    fid_wb = fopen(water_bottom,'r');
    wb_idx = fread(fid_wb,'float32');
    
    % zero pad data
    %sas_result_trc = [ava(1:ns,:);zeros(seismic.n_samples-ns,ntraces)];
    sas_result_trc_zero_pad = [sas_result_trc;zeros(seismic.n_samples-size(sas_result_trc,1),size(sas_result_trc,2))];
    
    %Unflatten data
for kk = 1:length(wb_idx)
    test(:,kk) = circshift(sas_result_trc_zero_pad(:,kk),wb_idx(kk));
end
    
    %ilxl_read = single(seismic.trace_ilxl_bytes(start_trace_idx:end_trace_idx,1:2));
    typecast(seismic.trace_ilxl_bytes(start_trace_idx:end_trace_idx,1:1),'int32');
    typecast(seismic.trace_ilxl_bytes(start_trace_idx:end_trace_idx,2:2),'int32');
    %ilxl_read = seismic.trace_ilxl_bytes(start_trace_idx:end_trace_idx,1:2);

    results_out{1,1} = 'ilxl numbers';
    %results_out{1,2} = {ilxl_read};
    results_out{1,2} = {seismic.trace_ilxl_bytes(start_trace_idx:end_trace_idx,1:2)};
    results_out{2,1} = 'sas_result_trc_ordered_unflattened';
    results_out{2,2} = test;

    sample_rate = 0.004;
    i_block = block;
    
    %sw_segy_write_clean(results_out,i_block, n_blocks, sample_rate, output_dir);
    sw_segy_write_lite_gathers(results_out,i_block, sample_rate, output_dir)
    fprintf('Block %d written to file \n',block)
end
   %% Remove temporary files created
   delete(strcat(output_dir,'sas_result_slice_*_trc_block_*.bin'));
    
    

end