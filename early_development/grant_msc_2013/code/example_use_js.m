% Example data read combining the matlab codes
% 'segy_make_structure_lite_tg.m', 'node_segy_read_traces_lite_tg.m' and
% 'trim_tg.m'. 
% the code begins by scanning the headers in the specified file and then
% extracting the desired information. The amplitude information is then
% read in using the code; 'node_segy_read_traces_lite_tg.m'

% Author: Troy Grant
% Date: 31/05/2013
%clear
%clc

%% Scan headers in file
filepath = '/apps/gsc/matlab-mcode-beta/grant_msc_2013/data/';
filename = 'gathers_starting_model_il1822_xl2100-3100_0-6000.segy';
%filename = 'starting_model_trim_sum_block_1.segy';
%for jj = 1:2 
il_byte = 189;
xl_byte = 193;
varargin = 37;
seismic = segy_make_structure_lite_tg(filepath,filename,il_byte,xl_byte,varargin);

%% Read the amplitude data
n_blocks = 1;
i_block = 1;
dec = 0;

[traces, ilxl_read] = node_segy_read_traces_lite_tg(seismic,i_block,n_blocks,dec);


%% Applying *whatever the code does* to the individual gathers within the matrix 
zsmooth = 40;
shift = 4;

count = 1;

n_gathers = length(unique(seismic.extra_bytes_data)); % works if all gathers same fold

results_out{2,2} = zeros(seismic.n_samples,seismic.n_traces/n_gathers);
results_out{3,2} = results_out{2,2};
results_out{4,2} = results_out{2,2};
for ii = 1:n_gathers:seismic.n_traces;
    data = traces(:,ii:(ii+39));
    [results_out{2,2}(:,count), results_out{3,2}(:,count), results_out{4,2}(:,count)] = trim_tg(data,shift,zsmooth);
    count = count + 1;
end

% Testing the array allocator for the trim-sum outputs - suggest altering
% so the row number in the array) changes with each application of a new file (e.g the
% starting model and then final model

results_out{1,1} = 'Some info here';
results_out{1,2} = % these should be inline crossline numbers
results_out{2,1} = 'trim_data';
results_out{3,1} = 'trim_shifts';
results_out{4,1} = 'trim_sum';

% Save result as segy
sw_segy_write_clean_tg(results_out,i_block,n_blocks,4,filepath,varargin);

% Run for both datasets
%filename = 'gathers_final_model_il1822_xl2100-3100_0-6000.segy';
%end 