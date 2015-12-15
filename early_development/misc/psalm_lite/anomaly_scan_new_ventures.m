% Scan segy
[seismic] = segy_make_structure_lite(filepath,filename,il_byte,xl_byte,varargin);

% Read segy
[traces, ilxl_read] = node_segy_read_traces_lite(seismic,i_block,n_blocks,dec);

% Track water bottom
[wb_idx] = water_bottom_flatten_lite(traces);

% Flatten volume
[traces_flat] = trace_flatten(traces,wb_idx);

% Sort to slice order


% Run anomaly scanner
seismic_anomaly_spotter_lite(seismic_mat_path,i_block,n_blocks,window_length,binary_slice_path,output_dir)

% Run connectivity

% Save results