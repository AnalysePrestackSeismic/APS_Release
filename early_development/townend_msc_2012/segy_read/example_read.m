segy = segy_read_files('/data/NOR/dtect/gullris_msc2012_townend/Misc');

% Read first 1000 traces
traces = segy_read_traces(segy{1},1,segy{1}.n_traces,0,0);

%figure(1)

%imagesc(traces.data)

