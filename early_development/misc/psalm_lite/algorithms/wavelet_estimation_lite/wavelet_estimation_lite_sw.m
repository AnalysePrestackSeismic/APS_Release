function [] = wavelet_estimation_lite_sw(seismic_mat_path,i_block,n_blocks,output_dir,seismic_filepath)

%%-------------------------------------------------------------------------
%   seismic_mat_path - path to flat binary file (.mat_lite) including
%   filename
%   i_block - set to '1' (in single quotes)
%   n_blocks - number of blocks to divide input into (also in
%   single quotes)
%   output_
%     seismic = load(seismic_mat_path);

%seismic.filepath = '/segy/KEN/2012_L10ab_final/filtered_angle_stacks/113j02_filtered_angle_stack_20-25_SEGY_disk.segy';
seismic.filepath = seismic_filepath;
%filename_index = regexp(seismic.filepath,'/','split');
%filename_index2 = regexp(filename_index(end),'(\d{2}-\d{2})','match');
%filename_index2 = regexp(seismic.filepath,'(\d{2}-\d{2})','match');
% filename_index = filename_index2{1};

tmp_seismic = fread(fopen(seismic_mat_path,'r'),'double');

%Collect number of samples, ILXL bytes and number of traces
seismic.n_samples = tmp_seismic(1,1);
seismic.trace_ilxl_bytes = reshape(tmp_seismic(2:end),3,[])';
seismic.n_traces = size(seismic.trace_ilxl_bytes,1);

%i_block = str2double(i_block);
n_blocks = str2double(n_blocks);

dec = 10;

%Read traces
[traces, ilxl_read] = node_segy_read_traces_lite(seismic,i_block,n_blocks,dec);

wb_idx = water_bottom_flatten_lite(traces);

n_traces_in_block = size(traces,2);

ns_win = 101;
ns_overlap = 50;
n_win = 1+floor((seismic.n_samples-ns_win)/(ns_win-ns_overlap));
w = zeros(2+ns_win,n_win);

win_idx = (0:1:ns_win-1)';

for ii = 1:n_win
    win_sub = bsxfun(@plus,wb_idx,win_idx+(ii-1)*ns_overlap);
    win_ind = bsxfun(@plus,win_sub,(0:seismic.n_samples:seismic.n_samples*(n_traces_in_block-1)));
    win_ind = win_ind(:,sum(win_sub > seismic.n_samples)==0);
    if isempty(win_ind)
        break
    end
    w(1,ii) = floor(ns_win/2)+(ii-1)*ns_overlap;
    w(2,ii) = size(win_ind,2);
    w(3:end,ii) = sum(abs(fft(traces(win_ind))),2);
end

save(strcat(char(output_dir),'wavelet_block_',num2str(i_block)),'w','-v7.3');
end

