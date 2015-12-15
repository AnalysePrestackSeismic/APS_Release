function [] = wavelet_estimation(job_meta_path,i_block)
%%-------------------------------------------------------------------------
% Estimate wavelets for DIGI
% Pick water bottom on 2/3 of max angle stack to get reliable pick
% Use this pick to flatten all other angle stacks
% Estimate wavelets from flattened angle stacks
% Input:
%
% Output:
%%
plot_on = 1;
% Load job meta information 
job_meta = load(job_meta_path);

% Wavelet estimation parameters
ns_win = 128;
ns_overlap = 32;

% Make directory to save results
if exist(strcat(job_meta.output_dir,'wavelets/'),'dir') == 0
    mkdir(strcat(job_meta.output_dir,'wavelets/'))
end
wav_directory = strcat(job_meta.output_dir,'wavelets/');

%%
% Read traces for 2/3 max angle stack
vol_index_wb = ceil(job_meta.nvols*0.6667);
[~, traces{vol_index_wb}, ilxl_read{vol_index_wb}] = ...
    node_segy_read(job_meta_path,num2str(vol_index_wb),i_block);

% Pick water bottom
[wb_idx] = water_bottom_picker(traces{vol_index_wb},0);
wb_idx(wb_idx < 0) = 1;
win_sub = bsxfun(@plus,wb_idx,(0:job_meta.n_samples{vol_index_wb}-max(wb_idx))');
win_ind = bsxfun(@plus,win_sub,(0:job_meta.n_samples{vol_index_wb}:...
job_meta.n_samples{vol_index_wb}*(size(traces{vol_index_wb},2)-1)));
%%
% Loop over all volumes and windows to estimate wavelets
for i_vol = 1:1:job_meta.nvols
    [~, traces, ~, ~] = node_segy_read(job_meta_path,num2str(i_vol),i_block);
    traces = traces(win_ind);
    [n_samples,n_traces] = size(traces);
 
    start_index = 1:ns_win-ns_overlap-1:n_samples-ns_win;
    end_index = start_index+ns_win-1;
    if end_index(end) < n_samples
        end_index(end+1) = n_samples;
        start_index(end+1) = n_samples-ns_win+1;
    end
    n_win = length(start_index);
    w = zeros(2+ns_win,n_win);
    for ii = 1:n_win        
        % Estimate wavelets and store meta information
        if ii == n_win
            w(1,ii) = end_index(ii);
        else
            w(1,ii) = start_index(ii);
        end
        w(2,ii) = n_traces;
        %NFFT = 2^nextpow2(size(win_sub,1));
        w(3:end,ii) = sum(abs(fft(traces...
            (start_index(ii):end_index(ii),:))),2); 
        % index as linear index, summed not averaged!
        
        % Produce plots showing the windows        
        if ii == 1 && plot_on == 1
            figure(1)
            imagesc(traces)
            hold on
            plot(repmat(w(1,ii),n_traces),'--')
            hold off
        elseif plot_on == 1
            hold on
            plot(repmat(w(1,ii),n_traces),'--')
            hold off
            %figure(2)
            %subplot(n_win,1,ii);
            %imagesc(traces(start_index(ii):end_index(ii),:));
        end
        
    end
    
    % Save estimated wavelets for this volume
    fid_wav = fopen(strcat(wav_directory,job_meta.volumes{i_vol},'_fft_wavelets_block_',i_block,'.bin'),'w');
    fwrite(fid_wav,w,'float32');
    fclose(fid_wav);
end
if str2double(i_block) == 1
    % Add processing information to job meta
    job_meta.wav_directory = wav_directory;
    job_meta.ns_win = ns_win;
    job_meta.ns_overlap = ns_overlap;
    job_meta.n_win = n_win;
    save(job_meta_path,'-struct','job_meta','-v7.3');
end
end