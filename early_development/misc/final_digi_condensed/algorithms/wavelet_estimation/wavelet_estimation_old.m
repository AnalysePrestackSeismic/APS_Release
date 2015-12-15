function [] = wavelet_estimation(job_meta_path,i_block)
%%-------------------------------------------------------------------------
% Estimate wavelets for DIGI
% Pick water bottom on 2/3 of max angle stack to get reliable pick
% Use this pick to flatten all other angle stacks
% Estimate wavelets from flattened angle stacks
% Input:
%
% Output:
plot_on = 1;
% Load job meta information 
job_meta = load(job_meta_path);

% Read traces for 2/3 max angle stack
vol_index_wb = ceil(job_meta.nvols*0.6667);
[~, traces{vol_index_wb}, ilxl_read{vol_index_wb}] = ...
    node_segy_read(job_meta_path,num2str(vol_index_wb),i_block);

% Pick water bottom
[wb_idx] = water_bottom_picker(traces{vol_index_wb},0);
wb_idx(wb_idx < 0) = 1;
% Load block for remaining angle stacks
for i_vol = 1:1:job_meta.nvols
    
    if i_vol ~= str2double(vol_index_wb) % don't repeat load for previously read stack
        % Read traces
        [~, traces{i_vol}, ilxl_read{i_vol}] = ...
        node_segy_read(job_meta_path,num2str(i_vol),i_block);
    end
end
%clearvars traces
% Estimate wavelets from traces block

% Check that all blocks have the same number of traces
for i_vol = 1:1:job_meta.nvols-1
    if size(traces{i_vol}) == size(traces{i_vol+1})
        n_samples = size(traces{1,i_vol},1);
        n_traces_in_block = size(traces{1,i_vol},2);
    else
        % Use ilxl_read to work out missing traces and correct for it
        k = 1;
    end
end

% Wavelet estimation parameters
ns_win = 128;
ns_overlap = 32;
n_win = 1+floor((n_samples-ns_win)/(ns_win-ns_overlap));
w = zeros(2+ns_win,n_win);
win_idx = (0:1:ns_win-1)';

% Make directory to save results
if exist(strcat(job_meta.output_dir,'wavelets/'),'dir') == 0
    mkdir(strcat(job_meta.output_dir,'wavelets/'))
end
wav_directory = strcat(job_meta.output_dir,'wavelets/');

% Loop over all volumes and windows to estimate wavelets
for i_vol = 1:1:job_meta.nvols
    for ii = 1:n_win
        win_sub = bsxfun(@plus,wb_idx,win_idx+(ii-1)*ns_overlap);
        % Make linear indices
        win_ind = bsxfun(@plus,win_sub,(0:n_samples:n_samples*(n_traces_in_block-1)));
        win_ind = win_ind(:,sum(win_sub > n_samples)==0);
        if isempty(win_ind)
            break
        end
        
        % Produce plots showing the windows        
        if ii == 1 && plot_on == 1
            figure(i_vol)
            imagesc(traces{i_vol})
            hold all
            %plot(win_sub(1,:))
            %plot(win_sub(end,:))
            plot(wb_idx+(floor(ns_win/2)+(ii-1)*ns_overlap),'--');
            hold off
        elseif plot_on == 1
            figure(i_vol)
            hold all
            %plot(win_sub(1,:))
            %plot(win_sub(end,:))
            plot(wb_idx+floor(ns_win/2)+(ii-1)*ns_overlap,'--');
            hold off
        end
        
        % Estimate wavelets and store meta information
        w(1,ii) = floor(ns_win/2)+(ii-1)*ns_overlap;
        w(2,ii) = size(win_ind,2);
        %NFFT = 2^nextpow2(size(win_sub,1));
        w(3:end,ii) = sum(abs(fft(traces{i_vol}(win_ind))),2); 
        % index as linear index, summed not averaged!
    end
    
    % Save estimated wavelets for this volume
    fid_wav = fopen(strcat(wav_directory,job_meta.volumes{i_vol},'_fft_wavelets_block_',i_block,'.bin'),'w');
    fwrite(fid_wav,w,'float32');
    w = zeros(2+ns_win,n_win);
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