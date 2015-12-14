function print_job_meta_info(job_meta_path)
%[seismic, traces, ilxl_read, offset_read] = node_segy_read(job_meta_path,vol_index,i_block)

job_meta = load(job_meta_path);

job_meta

end