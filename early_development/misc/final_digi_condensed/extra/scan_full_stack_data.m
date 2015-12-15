function scan_full_stack_data(filepath,filename,il_byte,xl_byte,output_location)

% Scan seismic file
seismic = segy_make_structure_lite(filepath,filename,il_byte,xl_byte);

n_blocks = 20;

% Distribute to cluster

for i_block = 1:1:n_blocks
    % Read volume
    [traces, ~, ~] = node_segy_read_traces_lite_gathers(seismic,i_block,n_blocks,dec);

    % Pick water bottom
    [wb_idx(:,i_block)] = water_bottom_flatten_lite(traces);
end

% Read or create a slice optimised version

% Run seismic anomaly spotter

% Save as segy

end