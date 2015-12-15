function [seismic, traces, ilxl_read, offset_read] = node_segy_read_traces_lite_gathers(seismic_mat_path,i_block,n_blocks,dec)

    % Load seismic structure (scanned trace headers)
    % seismic = load(seismic_mat_path);
    % Check with full stack or gathers
    
    tmp_seismic = fread(fopen(seismic_mat_path,'r'),'double');

    % read inline and crossline bytes from index file
    % read any other bytes too (e.g. vargargin for offset etc)
    seismic.filepath = char(tmp_seismic(1:2000,1)');
    seismic.file_type = tmp_seismic(2001,1);
    seismic.n_samples = tmp_seismic(2002,1);
    seismic.trace_ilxl_bytes = reshape(tmp_seismic(2003:end),3,[])';
    seismic.n_traces = size(seismic.trace_ilxl_bytes,1);
    
    
    % read the segy headers for use in output ############################
    % Read textual header (3200 bytes)
    seismic.text_header = fread(seismic.fid,3200,'uchar');
    seismic.text_header = char(ebcdic2ascii(reshape(seismic.text_header,80,40)'));
    
    % Read binary header as
    seismic.binary_header=fread(seismic.fid,400,'uint8');
    
    % Re-interpret binary header as uint16 or uint32 as required
    two_bytes=seismic.binary_header(1:2:399)*256+seismic.binary_header(2:2:400);
    four_bytes=((seismic.binary_header(1:4:9)*256+seismic.binary_header(2:4:10))*256+seismic.binary_header(3:4:11))*256+seismic.binary_header(4:4:12);
    seismic.binary_header=[four_bytes(1:3);two_bytes(7:200)];
    
    seismic.n_samples = seismic.binary_header(8);
    seismic.s_rate = seismic.binary_header(6);
    seismic.file_type = seismic.binary_header(10);
    
    %#####################################################################
    
    
        
    
    
    if seismic.file_type < 1 || seismic.file_type > 5
       seismic.file_type = input('No standard seismic file type. Please enter seismic file type (1 (IBM),2,3,4 or 5 (IEEE)): ', 's');
       seismic.file_type = str2num(seismic.file_type);
    end
    
    %%% PRE-STACK DATA %%% 
    % Scan first 1000 traces to check whether they are gathers   
    if length(unique(seismic.trace_ilxl_bytes(1:1000,1:2),'rows')) < 1000 || dec == -1;
        is_gathers = 1;
        % there are duplicate inline and crossline locations
        % likely to be a gather file so traces should output complete
        % gathers
        % unique(seismic.trace_ilxl_bytes(:,1:2),'rows')
        
        if i_block == 1
            start_trace_idx = 1+((i_block-1)*floor(seismic.n_traces/n_blocks));
        else
            start_trace_idx = 1+((i_block-1)*floor(seismic.n_traces/n_blocks));
            start_trace_idx_il_xl = seismic.trace_ilxl_bytes(start_trace_idx,1:2);
            
            while sum(start_trace_idx_il_xl == seismic.trace_ilxl_bytes(start_trace_idx,1:2)) == 2
                start_trace_idx = start_trace_idx + 1;   
            end
            % start_trace_idx = start_trace_idx + 1;
        end
        
        if i_block == n_blocks
            end_trace_idx = seismic.n_traces;
        else
            end_trace_idx = i_block*floor(seismic.n_traces/n_blocks);
            end_trace_idx_il_xl = seismic.trace_ilxl_bytes(end_trace_idx,1:2);

            while sum(end_trace_idx_il_xl == seismic.trace_ilxl_bytes(end_trace_idx,1:2)) == 2
                end_trace_idx = end_trace_idx + 1;   
            end
            end_trace_idx = end_trace_idx - 1;
        end  
        
    %%% STACK DATA %%%    
    else % likely to be a full stack volume
        % Calculate the trace to start reading from for this block
        is_gathers = 0;
        start_trace_idx = 1+((i_block-1)*floor(seismic.n_traces/n_blocks));

        % Calculate the last trace to read for this block (adjusts for extra traces on last block)
        if i_block == n_blocks
            end_trace_idx = seismic.n_traces;
        else
            end_trace_idx = i_block*floor(seismic.n_traces/n_blocks);
        end        
    end

    %%% READ THE DATA %%%
    % Calculate the number of traces to read in this block
    n_traces_to_read = end_trace_idx-start_trace_idx+1;

    % Open the seismic segy file
    fid = fopen(char(seismic.filepath),'r','b');
    
    fseek(fid,seismic.trace_ilxl_bytes(start_trace_idx,3)-240,'bof');

    if seismic.file_type == 1
    % Convert traces from IBM32FP read as UINT32 into IEEE64FP (doubles)
        %traces_tmp = fread(fid,[60+seismic.n_samples,n_traces_to_read],strcat(num2str(seismic.n_samples),'*uint32=>uint32'));
        traces_tmp = fread(fid,[60+seismic.n_samples,n_traces_to_read],'*uint32');
        ilxl_read = traces_tmp(48:49,:)'; % what happens if the inline and crossline are not in this location        
        offset_read = traces_tmp(10,:)';
        traces = (1-2*double(bitget(traces_tmp(61:end,:),32))).*16.^ ...
          (double(bitshift(bitand(traces_tmp(61:end,:),uint32(hex2dec('7f000000'))),-24))-64).* ...
          (double(bitand(traces_tmp(61:end,:),uint32(hex2dec('00ffffff'))))/2^24);
    elseif seismic.file_type == 2 
        disp('This seismic file type is not currently supported. Please speak to Charles Jones.');
    elseif seismic.file_type == 5
    % Traces are IEEE32FP (doubles)   
        traces = fread(fid,[60+seismic.n_samples,n_traces_to_read],strcat(num2str(seismic.n_samples),'*float32'));
        trace_headers = typecast(single(reshape(traces(1:60,:),1,60*n_traces_to_read)),'int32');  
        trace_headers = reshape(trace_headers,60,n_traces_to_read);
        ilxl_read = trace_headers(48:49,:)';        
        offset_read = trace_headers(10,:)';        
        traces = traces(61:end,:);            
    else
        disp('This seismic file type is not currently supported. Please speak to Charles Jones.');
    end

    fclose(fid);     
end

function [dec_ilxl_bytes] = decimate(seismic,dec)

    ils = unique(seismic.trace_ilxl_bytes(:,1));

    % Make a logical mask (e.g. 1;0;1;0;1;0...) to decimate inlines
    dec_il_mask = logical(repmat([1;zeros(dec-1,1)],ceil(length(ils)/dec),1));
    dec_il_mask = dec_il_mask(1:length(ils),1);
    
    % Decimate inlines
    dec_il = ils(dec_il_mask,1);
    
    % Find all xlines on the newly decimated inlines
    dec_ilxl_bytes = seismic.trace_ilxl_bytes(ismember(seismic.trace_ilxl_bytes(:,1),dec_il),:);
        
    % Make a xline number mask to decimate xlines
    dec_xl_mask = (min(seismic.trace_ilxl_bytes(:,2)):dec:max(seismic.trace_ilxl_bytes(:,2)))';

    % Decimate xlines
    dec_ilxl_bytes = dec_ilxl_bytes(ismember(dec_ilxl_bytes(:,2),dec_xl_mask),:);
end