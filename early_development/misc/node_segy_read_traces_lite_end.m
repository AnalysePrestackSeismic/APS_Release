function [] = node_segy_read_traces_lite_end()

for i_file = 1:1:14

    seismic.filepath = strcat('/data/CNS/segy/endeavour_gathers/file_j5456415616s',num2str(i_file));

    n_blocks = 5;

    for i_block = 1:1:n_blocks;
        write_data(seismic,i_block,n_blocks,0);    
    end
    
end

end

function [] = write_data(seismic,i_block,n_blocks,dec)

    seismic.n_samples = 1751;
    % Calculate number of traces    
    ll=dir(seismic.filepath);
    seismic.n_traces=0.25*(ll.bytes-3600)/(seismic.n_samples+60);   

    % Load seismic structre (scanned trace headers)
    % seismic = load(seismic_mat_path);
        
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
    
    %if i_block == 1;
        skip_textual_binary = 3600;
        trc_head = 240;
        trc_length = seismic.n_samples*4;
        seismic.trace_ilxl_bytes(1:1:seismic.n_traces,1) = ....
                     skip_textual_binary +...
                     cumsum(repmat(trc_head,1,seismic.n_traces)) +... 
                     (0:1:seismic.n_traces-1).*trc_length;
    %end
    
    % Open the seismic segy file
    fid = fopen(seismic.filepath,'r','l');

    % Skip to the first sample of the first trace to read for this block
    
    ebc = fread(fid,900,'float32');
    
    fclose(fid);
    % fseek(fid,3600,'bof');
    
    fid = fopen(seismic.filepath,'r','l');

    % Read all traces for this block
    fseek(fid,seismic.trace_ilxl_bytes(start_trace_idx,1)-240,'bof');   
    
    traces_tmp = fread(fid,[60+seismic.n_samples,n_traces_to_read],strcat(num2str(seismic.n_samples),'*uint32=>uint32','l'));
    ilxl_read = traces_tmp(48:49,:)';
    if dec ~= 0
        [trace_ilxl_bytes_tmp] = decimate(seismic,dec);
        dec_mask = ismember(ilxl_read,trace_ilxl_bytes_tmp(:,1:2),'rows')';
        traces_tmp = traces_tmp(:,dec_mask);
        ilxl_read = ilxl_read(dec_mask,:);
    end
    % Convert traces from IBM32FP read as UINT32 into IEEE64FP (doubles)
%     traces = (1-2*double(bitget(traces_tmp(61:end,:),32))).*16.^ ...
%       (double(bitshift(bitand(traces_tmp(61:end,:),uint32(hex2dec('7f000000'))),-24))-64).* ...
%       (double(bitand(traces_tmp(61:end,:),uint32(hex2dec('00ffffff'))))/2^24);
%            
%     fclose(fid);
    
    file_out1 = strcat(seismic.filepath,'_headers');
    fid = fopen(file_out1,'w','l');
    fwrite(fid,traces_tmp(1:60,:),'uint32','l');    
    fclose(fid);
    fid = fopen(file_out1,'r','b');
    traces_tmp(1:60,:) = fread(fid,[60,seismic.n_traces],'uint32','b');
    fclose(fid);
    
    tic
    file_out2 = strcat(seismic.filepath,'_block_',num2str(i_block),'_final');
    toc
    
    fid = fopen(file_out2,'w','l');
    fwrite(fid,ebc,'float32','l');    
    fclose(fid);
    fid = fopen(file_out2,'a','b');
    fwrite(fid,traces_tmp,'uint32','b');
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