% Example data read combining the matlab codes
% 'segy_make_structure_lite_tg.m', 'node_segy_read_traces_lite_tg.m',
% 'trim_tg.m' and 'sw_segy_write_clean_tg.m'.
% The code begins by scanning the headers in the specified file and then
% extracting the desired information. The amplitude information is then
% read in using the code; 'node_segy_read_traces_lite_tg.m' to convert
% everything from segy into a format that matlab can read (i.e a matrix).
% The gathers are then extracted and run through the file 'trim_tg.m' to
% calculate the amount of NMO. The output is then rewritten as a segy file.

% Author: Troy Grant
% Date: 31/05/2013
clear
clc

%% Scan headers in file
filepath = '/apps/gsc/matlab-mcode-beta/grant_msc_2013/data/';
filename = 'gathers_starting_model_il1822_xl2100-3100_0-6000.segy';
il_byte = 189;
xl_byte = 193;
varargin = 37;

for jj = 1:2;
seismic = segy_make_structure_lite_tg(filepath,filename,il_byte,xl_byte,varargin);
%% Read the amplitude data
n_blocks = 1;
i_block = 1;
dec = 0;
[traces, ilxl_read] = node_segy_read_traces_lite_tg(seismic,i_block,n_blocks,dec);

%% Calculating the amount of NMO between the individual gathers within the matrix 
% Defining the variables, counters and parameters needed 
zsmooth = 25;
shift = 5;
count = 1;
% trim_data = zeros(size(traces));
% trim_shift = zeros(size(traces));
fold = min(length(unique(seismic.extra_bytes_data))); % Smallest gather in the data
% trim_sum = zeros(seismic.n_samples,seismic.n_traces/fold); % Creating the trim_sum matrix full of zeros - by using the smallest
% fold in any of the gathers, we return the largest possible matrix of
% trim_sum with which to fill
no_traces = 0; % Number of traces within each gather
no_gathers = 0; % Number of gathers
prev_trace_no = 1; % Trace number in file 'traces' to begin the data calculation (changes with each subsequent gather)

% Creating the array results_out which is the input for the segy writing
% code. The array will consist of zeros that will be filled later
results_out{2,2} = zeros(seismic.n_samples,seismic.n_traces);
results_out{3,2} = results_out{2,2};
results_out{4,2} = zeros(seismic.n_samples,seismic.n_traces/fold);

% Calculating the fold of each gather and counting the number of gathers in
% the entire dataset
A = 0;
for ii = 2:length(seismic.extra_bytes_data);
    % A = seismic.extra_bytes_data(ii-1,:); % A in this instance is simply an arbitrary counter 
    B = seismic.extra_bytes_data(ii,:); % B in this instance is another arbitrary counter 
    no_traces = no_traces + 1;
    % If B < A this means we have reached the end of the gather - the
    % n_gathers parameter is updated and we run the gather through the code
    % 'trim_tg.m', writing everything to the corresponding files
    if B < A 
        no_gathers = no_gathers + 1;
        data = traces(:,prev_trace_no:no_traces);
        %[trim_data(:,prev_trace_no:no_traces), trim_shift(:,prev_trace_no:no_traces), trim_sum(:,count)] = trim_tg(data,shift,zsmooth);
        [results_out{2,2}(:,prev_trace_no:no_traces), results_out{3,2}(:,prev_trace_no:no_traces), results_out{4,2}(:,count)] = trim_tg(data,shift,zsmooth);
        ilxl(:,count) = seismic.trace_ilxl_bytes(ii-1,:); % Adding the inline, crossline and byte numbers for each gather in the file trim_sum
        % the format in descending row order is: inline, crossline, byte
        count = count + 1;
        prev_trace_no = no_traces + 1;
    end
    A = seismic.extra_bytes_data(ii,:); % A in this instance is simply an arbitrary counter 
end

%
results_out{4,2} = results_out{4,2}(:,1:no_gathers);

% Testing the array allocator for the trim-sum outputs - suggest altering
% so the row number in the array) changes with each application of a new file (e.g the
% starting model and then final model
results_out{1,1} = 'Inline and crossline numbers';
results_out{1,2}{1,1} =  ilxl_read;
results_out{1,2}{2,1} =  uint32(seismic.extra_bytes_data);
results_out{1,2}{1,2} =  ilxl_read;
results_out{1,2}{2,2} =  uint32(seismic.extra_bytes_data);
results_out{1,2}{1,3} =  uint32(ilxl(1:2,:)');
results_out{1,2}{2,3} =  uint32(ones(no_gathers,1));

if jj == 1
    results_out{2,1} = 'starting_model_trim_data';
    results_out{3,1} = 'starting_model_trim_shifts';
    results_out{4,1} = 'starting_model_trim_sum';
else
    results_out{2,1} = 'final_model_trim_data';
    results_out{3,1} = 'final_model_trim_shifts';
    results_out{4,1} = 'final_model_trim_sum';  
end

% Save result as segy
sw_segy_write_clean_tg(results_out,i_block,n_blocks,4,filepath,varargin);

% Run for both datasets
filename = 'gathers_final_model_il1822_xl2100-3100_0-6000.segy';
end 