% Example data read

% Scan headers in file
%filepath = '/media/30131145/msc/data/';
filepath = '/apps/gsc/matlab-mcode-beta/grant_msc_2013/data';
filename = 'gathers_final_model_il1822_xl2100-3100_0-6000.segy';
il_byte = 189;
xl_byte = 193;
varargin = 37;
seismic = segy_make_structure_lite_tg(filepath,filename,il_byte,xl_byte,varargin);

% Read the amplitude data
n_blocks = 1;
i_block = 1;
dec = 0;

[traces, ilxl_read] = node_segy_read_traces_lite_tg(seismic,i_block,n_blocks,dec);

% Do something with the data

%[trim_data, trim_shift] = trim_tg(data,shift,zsmooth);