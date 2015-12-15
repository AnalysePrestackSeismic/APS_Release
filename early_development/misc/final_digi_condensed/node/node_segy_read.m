function [seismic, traces, ilxl_read, offset_read] = node_segy_read(job_meta_path,vol_index,i_block)
%[seismic, traces, ilxl_read, offset_read] = node_segy_read(job_meta_path,vol_index,i_block)

i_block = str2double(i_block);
job_meta = load(job_meta_path);

if i_block > job_meta.n_blocks
    seismic = NaN;
    traces = NaN;
    ilxl_read = NaN;
    offset_read = NaN;
    fprintf('Job only has %d blocks...\n',job_meta.n_blocks)
    return
end

vol_keys = segy_index_byte_finder(job_meta_path,job_meta.block_keys(i_block,:),vol_index);

vol_index = str2double(vol_index);
vol_name = job_meta.volumes{vol_index};
    
% find the names of the blocks for the given volume
ii = 1;
for i_files = 1:1:size(job_meta.files,1)
    if strfind(job_meta.files{i_files},vol_name) == 1
        blocks{ii,1} = job_meta.files{i_files};
        ii = ii + 1;
    end
end    
   
% loop over the blocks
for ii_block = 1:1:size(blocks,1)
    if size(vol_keys{ii_block},1) > 2
        vol_keys{ii_block}(end+1,3) = 0;
        seismic = segy_read_binary(strcat(job_meta.paths{1},blocks{ii_block}));  

        bytes_per_sample = 4;
        trc_head = 240;
        trc_length = seismic.n_samples*bytes_per_sample; % cpuld add this to seismic structure

        s_key = 1;
        i_counter = 1;
        i_trace = 1;
        while s_key < size(vol_keys{ii_block},1)      
            if vol_keys{ii_block}(i_counter+1,3) - vol_keys{ii_block}(i_counter,3) == trc_length+240
                i_counter = i_counter + 1;
            else  % perform continuous read
                n_traces_to_read = i_counter-s_key+1;
                % Open the seismic segy file
                [traces{ii_block}(:,s_key:s_key+n_traces_to_read-1),...
                    ilxl_read{ii_block}(s_key:s_key+n_traces_to_read-1,:),...
                    offset_read{ii_block}(s_key:s_key+n_traces_to_read-1,:)] ...
                    = read_traces_segy(seismic,vol_keys{ii_block}(s_key,3)-240,n_traces_to_read);
                s_key = s_key + n_traces_to_read;
                i_counter = s_key;
            end
        end
        traces = cell2mat(traces);
        ilxl_read = cell2mat(ilxl_read);
        offset_read = cell2mat(offset_read);
    else
        % no traces
        seismic = segy_read_binary(strcat(job_meta.paths{1},blocks{ii_block}));
        traces = 0;
        ilxl_read = 0;
        offset_read = 0;
        % could to make zero entry
    end  
end
end
        
function [traces,ilxl_read,offset_read] = read_traces_segy(seismic,start_byte,n_traces_to_read)

    fid = fopen(char(seismic.filepath),'r','b');
    fseek(fid,start_byte,'bof');

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