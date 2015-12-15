function [] = wavelet_estimation_lite(seismic_mat_path,i_block,n_blocks,output_dir)

% declare variables ###################################################

i_block = str2double(i_block);
n_blocks = str2double(n_blocks);

% Wavelet estimation window parameters
ns_win = 101;
ns_overlap = 50;
gathyes = 10;

% Decimation factor
dec = 10;
n_blocks = (n_blocks*dec);

%filename_index = cell(1);
filename_index = 'testoutput3h';
%#######################################################################

% open the mat file and then get the segy file to read
%tmp_seismic = fread(fopen(seismic_mat_path,'r'),'double');

% read inline and crossline bytes from index file
% read any other bytes too (e.g. vargargin for offset etc)
%seismic.filepath = char(tmp_seismic(1:2000,1)');
%seismic.file_type = tmp_seismic(2001,1);
%seismic.n_samples = tmp_seismic(2002,1);
%seismic.trace_ilxl_bytes = reshape(tmp_seismic(2003:end),3,[])';
%seismic.n_traces = size(seismic.trace_ilxl_bytes,1);



% Wavelet Estimation section %
for i_block = 1:dec:n_blocks

    % Read traces from this block
    [seismic, traces, ilxl_read, ~] = node_segy_read_traces_lite_gathers(seismic_mat_path,i_block,n_blocks,gathyes);
     
    if strcmp('testoutput3h',filename_index) == 1;
        [filename_index] = regexp(seismic.filepath,'(\d{2}-\d{2})','match','once');
    end
    % Flatten to water bottom
    wb_idx = water_bottom_flatten_lite(traces);

    n_traces_in_block = size(traces,2);

    n_win = 1+floor((seismic.n_samples-ns_win)/(ns_win-ns_overlap));
    w = zeros(2+ns_win,n_win);

    win_idx = (0:1:ns_win-1)';

    % Wavelet estimation
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

    if i_block==1
        tmp_w = w;
    else
        tmp_w(2:end,:) = tmp_w(2:end,:)+w(2:end,:);
        if sum(logical(w(1,:))) > sum(logical(tmp_w(1,:)))
            tmp_w(1,:) = w(1,:);
        end
    end
        
end

% Average wavelets across blocks
avg_w = [tmp_w(1,:); bsxfun(@rdivide,tmp_w(3:end,:),tmp_w(2,:))];
avg_w = avg_w(:,logical(1-logical(sum(isnan(avg_w)))));

avg_w = avg_w(:,logical(sum(avg_w(2:end,:))));
avg_w_time = [avg_w(1,:); circshift(ifft(avg_w(2:end,:),'symmetric'),floor(ns_win/2))];

% Save average wavelets
if strcmp('',filename_index) == 1;
    filename_index = 'unknown_angle_range';
end    
save(strcat(char(output_dir),'average_wavelets_freq_',filename_index),'avg_w','-v7.3');
save(strcat(char(output_dir),'average_wavelets_time_',filename_index),'avg_w_time','-v7.3');
%save(strcat(char(output_dir),'average_wavelets_freq_',filename_index{1}{1}),'avg_w','-v7.3');
%save(strcat(char(output_dir),'average_wavelets_time_',filename_index{1}{1}),'avg_w_time','-v7.3');

% Compile final wavelets into cell array to be used in IG inversion
%all_wavelets_freq{1,a} = avg_w;
%all_wavelets_time{1,a} = avg_w_time;


% Save final wavelet files
%save(strcat(char(output_dir),'all_wavelets_freq'),'all_wavelets_freq','-v7.3');
%save(strcat(char(output_dir),'all_wavelets_time'),'all_wavelets_time','-v7.3');
end

% Remove temporary files created
% delete(strcat(output_dir,'wavelet_block_*.mat'));

% Plot Wavelets
%figure(1)

%for a = 1: size(angle_stacks,1)
%    
%     if a == 1
%         load('/data/Global/dtect/SRW_test_dataset/SEGY/matlab/test/average_wavelet_time_05-10.mat');
%         subplot(8,1,1);
%         imagesc(bsxfun(@rdivide,avg_w_time(2:end,:),max(avg_w_time(2:end,:))));
%     else
%         load(sprintf('/data/Global/dtect/SRW_test_dataset/SEGY/matlab/test/average_wavelet_time_%d-%d',a*5,a*5+5));
%         subplot(8,1,a);
%         imagesc(bsxfun(@rdivide,avg_w_time(2:end,:),max(avg_w_time(2:end,:))));
%     end
%     
%     
% end
