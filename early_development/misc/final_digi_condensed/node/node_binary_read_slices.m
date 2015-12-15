function [slices] = node_binary_read_slices(job_meta_path,...
    binary_slice_path,i_vol,i_slab,n_slab,window_length,varargin)

job_meta = load(job_meta_path);

% Count the number of slice ordered binary files in the path binary_slice_path
binary_slice_path = strcat(job_meta.paths{str2double(i_vol)},binary_slice_path);
[files_in,n_files] = directory_scan({binary_slice_path},'bin');
files_in.names = sort_nat(files_in.names);
% Get ns
file_to_open = strcat(binary_slice_path,files_in.names{1});
fid = fopen(file_to_open,'r');    
n_samples = fread(fid,1,'float32');
fclose(fid);

% Calculate the slice to start reading from for this slab

start_slice_idx = (1+((str2double(i_slab)-1)*floor(n_samples/str2double(n_slab))))-str2double(window_length);
%   start_slice_idx = (1+((i_block-1)*floor(ns/n_blocks)))-window_length;

% Calculate the last trace to read for this slab (adjusts for extra slices in last slab)
if str2double(i_slab) == str2double(n_slab)
    end_slice_idx = window_length+n_samples;
    %     end_slice_idx = window_length+ns;
else
    end_slice_idx = str2double(window_length)+(str2double(i_slab)*floor(n_samples/str2double(n_slab)));
    %     end_slice_idx = window_length+(i_block*floor(ns/n_blocks));
    
end

pad_start = [];
pad_end = [];

% Correct for window overshoots
if start_slice_idx <= 0
    pad_start = zeros(job_meta.vol_traces,1+abs(start_slice_idx));
    %pad_start = zeros(seismic.n_traces/n_blocks,1+abs(start_slice_idx));
    start_slice_idx = 1;
end
if end_slice_idx >= n_samples
    % if end_slice_idx >= ns
    pad_end = zeros(job_meta.vol_traces,end_slice_idx-n_samples);
    end_slice_idx = n_samples;
    %     end_slice_idx = ns;
end

% Calculate the number of slices to read in this slab
n_slices_to_read = end_slice_idx-start_slice_idx+1;

% Loop through all binary files and read the slices in this slab
% n_traces = 0;
for i_file = 1:n_files
    % Open the binary file
    file_to_open = strcat(binary_slice_path,files_in.names{i_file});
    fid = fopen(file_to_open,'r');
    
    n_samples = fread(fid,1,'float32');
    % ns = 600;
    
    % Calculate the number of traces in this binary file
    ll = dir(file_to_open);
    n_traces_in_file = (ll.bytes-4)/(4*n_samples);
    % n_traces(ii) = n_traces_in_file;
    
    % Skip to the first sample of the first slice to read for this slab
    fseek(fid,4+(4*n_traces_in_file*(start_slice_idx-1)),'bof');
    % Read all samples for this slab
    slices{i_file,1} = fread(fid,[n_traces_in_file,n_slices_to_read],'float32');
    
%     test_array(ii,1) = n_traces_in_file;
%     test_array(ii,2) = n_slices_to_read;
    
    fclose(fid);
    % fprintf('Slab read %d \r\n', ii);
end

% Convert from cell array to matrix (each column is a time slice)
slices = [pad_start,cell2mat(slices),pad_end];

if size(varargin) > 0;
    for ii = 1:1:size(slices(:,2)); slice = slices(:,ii); slice = reshape(slice,length(unique(seismic.trace_ilxl_bytes(:,2))),length(unique(seismic.trace_ilxl_bytes(:,1)))); imagesc(slice); title(num2str(ii)); pause(1); end
end

%%Section to test size of slices - compare to available memory.

% size of slices variable in MB
% slice_size = (getfield(whos('slices'),'bytes'))/1024^2

end



%% TO DO:

% NS is written into the flat binary files in the inversion code but is not
% yet being picked up properly here. Need to change this and make the file
% rely on NS rather than seismic.n_samples. Check pad_start is actually the
% right dimensions.