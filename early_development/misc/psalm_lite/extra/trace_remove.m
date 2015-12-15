function [] = trace_remove(seismic_filepath,seismic_mat_path,i_block,n_blocks,ilxl_remove)
% Remove traces specified by the list of il/xl locations 'ilxl_remove' from
% the segy rev1 file 'seismic_filepath'. The segy must first be scaned
% using segy_make_structure_lite (or similar) to produce the segy lookup
% table 'seismic_mat_path'. For files larger than the available memory,
% loop the function over 'n_blocks', incrementing the loop by 'i_block'.

    % Open and read the scanned segy file
    fid_scan = fopen(seismic_mat_path,'r');
    tmp_seismic = fread(fid_scan,'double');
    
    % Get the required info about the segy file
    seismic.n_samples = tmp_seismic(1,1);
    seismic.trace_ilxl_bytes = reshape(tmp_seismic(2:end),3,[])';
    seismic.n_traces = size(seismic.trace_ilxl_bytes,1);
    
    if i_block == 1
        n_traces_remove = sum(ismember(seismic.trace_ilxl_bytes(:,1:2),ilxl_remove,'rows'));
        fprintf('\nNumber to traces that will be removed is %d\n',n_traces_remove);
        ll = dir(seismic_filepath);
        fprintf('Output filesize will be %ld bytes\n',ll.bytes-n_traces_remove*(240+seismic.n_samples*4));
        clearvars n_traces_remove ll
    end
    
    % Tidy up
    clearvars tmp_seismic;
    fclose(fid_scan);
        
    % Calculate the trace to start reading from for this block
    start_trace_idx = 1+((i_block-1)*floor(seismic.n_traces/n_blocks));

    % Calculate the last trace to read for this block (adjusts for extra traces on last block)
    if i_block == n_blocks
        end_trace_idx = seismic.n_traces;
    else
        end_trace_idx = i_block*floor(seismic.n_traces/n_blocks);
    end

    % Calculate the number of traces to read in this block
    n_traces_to_read = end_trace_idx-start_trace_idx+1;

    % Open the seismic segy file
    fid_in = fopen(seismic_filepath,'r','b');

    % Read file header if on the first block
    if i_block == 1
        header = fread(fid_in,[900,1],'*uint32','b');   
    end
    
    % Skip to the first sample of the first trace to read for this block
    fseek(fid_in,seismic.trace_ilxl_bytes(start_trace_idx,3)-240,'bof');
    
    % Read all traces for this block
    fprintf('\nReading the %d traces which make up block %d of %d\n',n_traces_to_read,i_block,n_blocks);
    traces = fread(fid_in,[60+seismic.n_samples,n_traces_to_read],'*uint32','b');
    fclose(fid_in);
    
    % Rows 48 and 49 of traces are the ils/xls for segy rev1
    ilxl_read = traces(48:49,:)';

    % Make logical index for traces that need removing
    fprintf('Checking if any traces need removing... ')
    ilxl_remove_idx = ismember(ilxl_read,ilxl_remove,'rows');
    fprintf('found %d to remove\n',sum(ilxl_remove_idx));
    
    % Open/create the output file
    fid_out = fopen(strcat(seismic_filepath,'_trimmed'),'a','b');
    
    % Write the file header if on the first block
    if i_block == 1 
        fwrite(fid_out,header,'uint32','b');
    end
    
    % Write the traces headers and trace data for all traces except those in the remove list
    fprintf('Writing the %d traces which make up block %d of %d to file %s\n',...
        n_traces_to_read-sum(ilxl_remove_idx),i_block,n_blocks,strcat(seismic_filepath,'_trimmed'));
    if sum(ilxl_remove_idx) ~= 0
        traces = traces(:,~ilxl_remove_idx);
    end
    fwrite(fid_out,traces,'uint32','b');
    
    % Tidy up
    fclose(fid_out);
end


