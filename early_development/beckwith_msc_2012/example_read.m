segy = segy_read_files('/data/TZA/dtect/bg_tza_site_survey_vol2_msc/');

cutoff=0

traces = segy_read_traces(segy{1},1,segy{1}.n_traces,0,0,cutoff);

