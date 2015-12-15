function seismic = segy_read_binary(seismic_mat_path)
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
    
    fid = fopen(seismic_mat_path,'r');
    message = ferror(fid); 
    tmp_seismic = fread(fid,'double');
    seismic.filepath = char(tmp_seismic(1:filepath_length,1)');
    seismic.file_type = tmp_seismic(file_type_length,1);
    seismic.s_rate = tmp_seismic(s_rate_length,1);
    seismic.n_samples = tmp_seismic(n_samples_length,1);
    seismic.n_traces = tmp_seismic(n_traces_length,1);
    seismic.pkey = tmp_seismic(pkey_length,1);
    seismic.skey = tmp_seismic(skey_length,1);
    seismic.tkey = tmp_seismic(tkey_length,1);
    
    seismic.is_gather = tmp_seismic(is_gather_length,1);

    if seismic.is_gather == 1
        seismic.trace_ilxl_bytes = reshape(tmp_seismic(traces_length:end),8,[])';
    elseif seismic.is_gather == 0
        seismic.trace_ilxl_bytes = reshape(tmp_seismic(traces_length:end),5,[])';
    else

    end

    fclose(fid);
end