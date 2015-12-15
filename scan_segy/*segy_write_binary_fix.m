function seismic = segy_write_binary_fix(seismic_mat_path, file_mat)
%% Defination: 
% Input:
% Output:
% Writes to Disk:
%%
    % Seismic binary structure
    filepath_length = 2000;
    file_type_length = filepath_length+1;
    s_rate_length = file_type_length+1;
    n_samples_length = s_rate_length+1;
    n_traces_length = n_samples_length+1; 
    pkey_length = n_traces_length+1;
    skey_length = pkey_length+1;
    tkey_length = skey_length+1;
    is_gather_length = tkey_length+1;
    traces_length = is_gather_length+1;
    %%    
    fid = fopen(seismic_mat_path,'r');
    %message = ferror(fid);     
    fid_write = fopen(file_mat,'w'); 
    %message = ferror(fid_write);        
    tmp_seismic = fread(fid,'double');
    seismic.filepath = char(tmp_seismic(1:filepath_length,1)');
    seismic.filepath = deblank(seismic.filepath);
    
    %%
    % forgot to put /data infront of a path   
    seismic.filepath = regexprep(seismic.filepath, '/URY', '/data/URY');
    
    filepath_binary = uint64(seismic.filepath);
    pad_filepath = zeros(1,(2000-length(filepath_binary)));
    filepath_binary = [filepath_binary,pad_filepath];
     
    
    seismic.file_type = tmp_seismic(file_type_length,1);
    seismic.s_rate = tmp_seismic(s_rate_length,1);
    seismic.n_samples = tmp_seismic(n_samples_length,1);
    seismic.n_traces = tmp_seismic(n_traces_length,1);
    seismic.pkey = tmp_seismic(pkey_length,1);
    seismic.skey = tmp_seismic(skey_length,1);
    seismic.tkey = tmp_seismic(tkey_length,1);
    seismic.is_gather = tmp_seismic(is_gather_length,1);

    
% do not want to manipulate the binary structure    
%     if seismic.is_gather == 1
%         seismic.trace_ilxl_bytes = reshape(tmp_seismic(traces_length:end),8,[])';
%     elseif seismic.is_gather == 0
%         seismic.trace_ilxl_bytes = reshape(tmp_seismic(traces_length:end),5,[])';
%     else
% 
%     end
    seismic.trace_ilxl_bytes = tmp_seismic(traces_length:end);
    
    %%
    fwrite(fid_write,[filepath_binary';seismic.file_type;seismic.s_rate;seismic.n_samples;seismic.n_traces;seismic.pkey;seismic.skey;seismic.tkey;seismic.is_gather;seismic.trace_ilxl_bytes],'double');
    

    fclose('all');
end