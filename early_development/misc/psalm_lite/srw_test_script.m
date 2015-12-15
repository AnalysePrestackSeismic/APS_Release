seismic_mat_path = '/data/Global/dtect/SRW_test_dataset/SEGY/matlab/full_run/20-25_angle_stack_SRW.mat_lite';
binary_slice_path = '/data/Global/dtect/SRW_test_dataset/SEGY/matlab/full_run/digi_results/';

tmp_seismic = fread(fopen(seismic_mat_path,'r'),'double');
seismic.n_samples = tmp_seismic(1,1);
seismic.trace_ilxl_bytes = reshape(tmp_seismic(2:end),3,[])';
seismic.n_traces = size(seismic.trace_ilxl_bytes,1);
i_block = 1;
n_blocks = 100;
window_length = 50;

file_to_open = sprintf('%sdigi_minimum_energy_eer_projection_slices_block_%d.bin',binary_slice_path,i_block);
fid = fopen(file_to_open,'r');    
ns = fread(fid,1,'float32');
fclose(fid);

% Count the number of slice ordered binary files in the path binary_slice_path
count_files = sprintf('ls -B %sdigi_minimum_energy_eer_projection_slices_block_*.bin | wc -l',binary_slice_path);
[~,n_files] = system(count_files);
%n_files = str2double(n_files);
n_files = 100;

for seismic_n_traces = 1000:1000:10000
    count = seismic_n_traces-999;
    
start_slice_idx = (1+((i_block-1)*floor(ns/n_blocks)))-window_length;
%   start_slice_idx = (1+((i_block-1)*floor(ns/n_blocks)))-window_length;

% Calculate the last trace to read for this slab (adjusts for extra slices in last slab)
if i_block == n_blocks
    end_slice_idx = window_length+ns;
    %     end_slice_idx = window_length+ns;
else
    end_slice_idx = window_length+(i_block*floor(ns/n_blocks));
    %     end_slice_idx = window_length+(i_block*floor(ns/n_blocks));
    
end

pad_start = [];
pad_end = [];

% Correct for window overshoots
if start_slice_idx <= 0
 pad_start = zeros(seismic_n_traces,1+abs(start_slice_idx));
 start_slice_idx = 1;
end
if end_slice_idx >= ns
    % if end_slice_idx >= ns
pad_end = zeros(seismic_n_traces,end_slice_idx-ns);
end_slice_idx = ns;
    %     end_slice_idx = ns;
end

% Calculate the number of slices to read in this slab
n_slices_to_read = end_slice_idx-start_slice_idx+1;

% Loop through all binary files and read the slices in this slab
% n_traces = 0;
for ii = 1:100
    % Open the binary file
    file_to_open = sprintf('%sdigi_minimum_energy_eer_projection_slices_block_%d.bin',binary_slice_path,ii);
    fid = fopen(file_to_open,'r');
    
    ns = 600;
    % ns = 600;
    
    % Calculate the number of traces in this binary file
    %ll = dir(file_to_open);i_block
    %n_traces_in_file = (ll.bytes-4)/(4*ns);
    n_traces_in_file = floor(seismic_n_traces/n_blocks);
    %n_traces(ii) = n_traces_in_file;
    
    % Skip to the first sample of the first slice to read for this slab
    fseek(fid,4+(4*n_traces_in_file*(start_slice_idx-1)),'bof');
    % Read all samples for this slab
    slices{ii,1} = fread(fid,[n_traces_in_file,n_slices_to_read],'float32');
    
%     test_array(ii,1) = n_traces_in_file;
%     test_array(ii,2) = n_slices_to_read;
    
    fclose(fid);
    % fprintf('Slab read %d \r\n', ii);
end

% Convert from cell array to matrix (each column is a time slice)
%slices = [pad_start,cell2mat(slices),pad_end];
slices = cell2mat(slices);

%% - - - - - - - - - - - - - - -- - - - - - - - --- - - -- - - - -- -  - -


% Get some information about the seimsic data
tmp_seismic = fread(fopen(seismic_mat_path,'r'),'double');
seismic.n_samples = tmp_seismic(1,1);
seismic.trace_ilxl_bytes = reshape(tmp_seismic(2:end),3,[])';
seismic.n_traces = size(seismic.trace_ilxl_bytes,1);

% Work out the min/max ils/xls for the survey bounding box
ilmin = min(seismic.trace_ilxl_bytes(:,1));
ilmax = max(seismic.trace_ilxl_bytes(:,1));
ilinc = mode(diff(unique(seismic.trace_ilxl_bytes(:,1))));
iln = 1+((ilmax-ilmin)/ilinc);

xlmin = min(seismic.trace_ilxl_bytes(:,2));
xlmax = max(seismic.trace_ilxl_bytes(:,2));
xlinc = mode(diff(unique(seismic.trace_ilxl_bytes(:,2))));
xln = 1+((xlmax-xlmin)/xlinc);

% Make ils/xls of the survey bounding box
bounding_box = [reshape(repmat((ilmin:ilinc:ilmax),xln,1),[],1),repmat((xlmin:xlinc:xlmax)',iln,1),(1:1:iln*xln)'];

% Get the linear indices of the lives traces in the survey bounding box
live_traces = bounding_box(ismember(bounding_box(:,1:2),seismic.trace_ilxl_bytes(:,1:2),'rows'),3);






% Work out the number of windows to loop through
n_windows = size(slices,2)-2*window_length;

% Pre-allocate memory for the output zero-padded data probability slab
dataprob = zeros(iln*xln,n_windows);

for ii = 1:n_windows
    % Calculate the cdf
    data = slices(:,ii:ii+2*window_length);
    
    %------------------------------------------------
    % Temp section to store data size
    
    data_size(ii,1) = ii;
    data_size(ii,2) = (getfield(whos('data'),'bytes'))/1024^2;
    
    
    data = data(:);
    if min(data) ~= max(data)
        data = data(data ~= 0);                                 % remove hard zeros
        data = data(~isnan(data));                              % remove NaN's
        numdata = length(data);                                 % total number of real datapoints (without zeros and NaN's)
        widthbin = 3.5*std(data)/(numdata^(1/3));               % gives optimum bin width for gaussian distribution (fairly robust for other distributions)
        bins = ((min(data):widthbin:max(data)))';               % defines bin centres based on bin width and data extremes
        data_tmp = interp1(bins,bins,data,'nearest','extrap');  % interpolates data onto closest bin center
        check_bins = ~ismember(bins,data_tmp);                  % finds empty bins
        data_tmp = sort([data_tmp;bins(check_bins)]);           % temporarily pads data with empty bins to ensure no bins are empty, and sorts low to high
        [~,datapdf_tmp,~] = unique(data_tmp);                   % returns the highest index of the unique values in data
        datapdf_tmp = diff([0;datapdf_tmp]);                    % difference between indices is the number of values in a bin
        datapdf_tmp(check_bins) = 0;                            % reset the padded bins to zero
        datapdf = datapdf_tmp;                                  % put into a cell array
        datacdf = cumsum(datapdf)-datapdf/2;                    % cumulative sum across bins and re-centers to bin center
        datacdf = datacdf/max(datacdf);                         % normalises max of CDF to 1
        
        % Calculate probabilities from cdf
        data = slices(:,ii+window_length);                      % get central slice in window
        data = interp1(bins,bins,data,'nearest','extrap');      % interpolates data onto closest bin center
        [~,idxdataprob] = ismember(data,bins);                  % finds the bin index that the data is in
        dataprob(live_traces,ii) = 1-datacdf(idxdataprob);      % lookup the probability scores from the cdf
    end
end
















seismic_n_traces
variable_size = struct2cell(whos);
variable_size = (sum(cell2mat(variable_size(3,:))))/1024^2

slices = {0};
end