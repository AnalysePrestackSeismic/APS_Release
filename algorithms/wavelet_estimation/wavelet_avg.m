function [] = wavelet_avg(job_meta_path)
%% ------------------ Disclaimer  ------------------
% 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) makes no representation or warranty, express or implied, in 
% respect to the quality, accuracy or usefulness of this repository. The code
% is this repository is supplied with the explicit understanding and 
% agreement of recipient that any action taken or expenditure made by 
% recipient based on its examination, evaluation, interpretation or use is 
% at its own risk and responsibility.
% 
% No representation or warranty, express or implied, is or will be made in 
% relation to the accuracy or completeness of the information in this 
% repository and no responsibility or liability is or will be accepted by 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) in relation to it.

%% ------------------ License  ------------------ 
% GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
%% github
% https://github.com/AnalysePrestackSeismic/
%% ------------------ FUNCTION DEFINITION ---------------------------------
% wavelet_avg: function to create wavelet set for use in DIGI. Currently
% wavelets vary with time and angle. Spatially variant wavelets not yet
% implemented.
%   Arguments:
%	job_meta_path = path to .mat file created using
%       segy_make_job function.
%
%   Outputs:
%       all_wavelets_time.mat = mat file to be input into DIGI.
%       This .mat file contains .... ??
%       all_wavelets_freq (...), all_wavelets_time (...), max_wavelet_zgrid (...),
%       min_wavelet_zgrid(...)
%
%   Writes to Disk:
%       nothing 

job_meta = load(job_meta_path);             % Load job meta information 

if job_meta.is_gather == 0                  % For the case of using angle stacks
    i_vol_max = job_meta.nvols;             % Number of angle stacks
else                                        % For the case of using angle gathers
    i_vol_max = ((job_meta.tkey_max -  job_meta.tkey_min )/ job_meta.tkey_inc) +  1; % Number of angles in angle gather
end

%% Wavelet Database Management for all bocks-- preconditioning for averaging

loopfin = size(job_meta.liveblocks,1);      % Number of Live Blocks
lpi = 1;                                    % Loop Index
count = 1;                                  % Counter for blocks that have a corresponding wavelet file
while lpi <= loopfin
    i_block = job_meta.liveblocks(lpi);     % Reference the block numbers for live blocks
    tmpfileout = sprintf('opening file %s',strcat(job_meta.wav_directory,'fft_wavelets_block_',num2str(i_block),'.bin'));
  
    if exist(strcat(job_meta.wav_directory,'fft_wavelets_block_',num2str(i_block),'.bin'),'file') % check if the file exists
        disp(tmpfileout)
        fid_wav = fopen(strcat(job_meta.wav_directory,...
            'fft_wavelets_block_',num2str(i_block),'.bin'));    % open the binary wavelet file for the block
        w = fread(fid_wav,'float32');                           % Read the binary file and convert to 32 bit floating point precission
        n_vol = w(1);                                           % Number of volumes (angle stacks)
        n_win = w(2);                                           % Number of windows for wavelet estimation
        stdev = w(3);                                           % Standard Deviation :might like to think about calculating this properly based on live samples
        live_offset = w(4);                                     % The live angle volume
        wb_z_avg = w(5);                                        % Average Z waterbottom in the block?
        wave_offset = 6;                                        % This is the index in the wavlet file from where the individual wavelet data starts
        job_meta.stdev(i_block,1) = stdev;                      % Store in job meta file for current live block  
        job_meta.wb_z_avg(i_block,1) = wb_z_avg;                % Store in job meta file for current live block 
        job_meta.live_offset(i_block,1) = live_offset;          % Store in job meta file for current live block         
        fclose(fid_wav);                                        % Close file

        if count == 1                                           % Initialize the wavelet matrix if this is the first live block
            w_all = zeros(n_win*(job_meta.ns_win+2),n_vol,loopfin); % Matrix for storing all wavelets for all volumes for all live blocks, (Note: length of wavlet array = number of samples in wavelet) + 2
            w_all(:,:,count) = reshape(w(wave_offset:end),[],n_vol);% Reshape the Matrix sperating the different angle volumes in different columns for the first block              
            max_n_win = n_win;                                  % Intitialize maximum number of wavelet windows in block as number of windows in first block
            max_n_vol = n_vol;                                  % Intitialize maximum number of angle volumes in block as number of windows in first block
        else                                                    % If this is not the first live block, append already initialized wavelet matrix
            if n_win > max_n_win                                % If you encounter a block with bigger window than what has been encountered ljust append a slab for the extra window
                w_append = zeros((n_win-max_n_win)*(job_meta.ns_win+2),n_vol,loopfin); % Create new slab of the extra length of window for all volumes for all live blocks
                w_all = [w_all; w_append];                      % Append the new slab for sccomodating the extra length of window
                w_all(1:n_win*(job_meta.ns_win+2),:,count) = reshape(w(wave_offset:end),[],n_vol); % Reshape the Matrix sperating the different angle volumes in different columns for the current block 
                max_n_win = n_win;
            else                                                % If this block has the same or less number of wavelet windows than encountered before 
                w_all(1:n_win*(job_meta.ns_win+2),:,count) = reshape(w(wave_offset:end),[],n_vol);% Reshape the Matrix sperating the different angle volumes in different columns for the current block
            end
        end
        lpi = lpi + 1;                                          % Increment loop index 
        count  = count + 1;                                     % Increment counter
    else                                                        % If file doesnot exist
        tmpfileout = sprintf('Error opening file %s',strcat(job_meta.wav_directory,'fft_wavelets_block_',num2str(i_block),'.bin'));
        disp(tmpfileout)                                        % Display error in opening file
        lpi = lpi + 1;                                          % Increment loop index 
    end
end

% Note: Description of w_all:
% Dimension 1 - time windows of wavelets
% Dimension 2 - different volumes or angles (offsets)
% Dimension 3 - different blocks (live ones)
                    
%% Wavelet averaging  
    %need to take account of zeros
    w_sum = sum(w_all,3);                                       % Sum the wavelet spectrum including number of live traces in the blocks for all the blocks
    get_windows = squeeze(w_all(1:job_meta.ns_win+2:end,1,:));  % Create a matrix of window starts for all blocks (1 column for each block) 
    get_windows = unique(get_windows(get_windows ~= 0));        % Create an array of window starts till the maximum Z encountered in any block
    w_sum(1:job_meta.ns_win+2:size(w_sum,1),:) = repmat(get_windows,1,i_vol_max);
    w_sum = reshape(w_sum,job_meta.ns_win+2,[]);
    
    %max_traces = max(w_sum(2,:));
    %w_sum = w_sum(:,w_sum(2,:) == max_traces);
    
    w_avg = [w_sum(2,:); bsxfun(@rdivide,w_sum(3:end,:),w_sum(2,:))];                                   % Divide the summation of spectrums by total number of live traces  for all windows
    w_avg_time = [w_sum(1,:); circshift(ifft(w_avg(2:end,:),'symmetric'),floor(job_meta.ns_win/2))];    % Inverse fourier transform to get the time waveletl
    w_avg_freq = [w_sum(1,:); w_avg(2:end,:)];                                                          % convert wavelets from frequency to time
    
% perform averaging over smaller blocks    
  
%         ilxl_grid = job_meta.block_keys(job_meta.liveblocks,:)
%         ilxl_aperture_step = 2;
%         n_pos = (2*ilxl_aperture+1)^2;  
%         positions_right = 0:ilxl_aperture_step:ilxl_aperture;
%         positions_left = fliplr(-ilxl_aperture_step:-ilxl_aperture_step:-ilxl_aperture);    
%         pos_diff_xl = repmat([positions_left positions_right],sqrt(n_pos),1);
%         pos_diff_il = pos_diff_xl';
%         pos_diff_il =  pos_diff_il(:).*mode(job_meta.pkey_inc); % turn into vector
%         pos_diff_xl =  pos_diff_xl(:).*mode(job_meta.skey_inc);
%         rep_il_pos = bsxfun(@plus,ilxl_grid(:,1), ...
%             pos_diff_il');            
%         rep_xl_pos = bsxfun(@plus,ilxl_grid(:,2), ...
%             pos_diff_xl');
%         rep_il_pos = rep_il_pos';
%         rep_xl_pos = rep_xl_pos';
%         rep_il_pos = rep_il_pos(:);
%         rep_xl_pos = rep_xl_pos(:);
%         temp_positions = [rep_il_pos rep_xl_pos];
%     
% 
%     % Average wavelets across blocks
%     avg_w = [tmp_w(1,:); bsxfun(@rdivide,tmp_w(3:end,:),tmp_w(2,:))];
%     avg_w = avg_w(:,logical(1-logical(sum(isnan(avg_w)))));
% 
%     avg_w = avg_w(:,logical(sum(avg_w(2:end,:))));
    
    % Save average wavelets
    
    %save(strcat(job_meta.wav_directory,job_meta.volumes{i_vol},'_avg_w_freq.mat'),'avg_w','-v7.3');
    %save(strcat(job_meta.wav_directory,job_meta.volumes{i_vol},'_avg_w_freq.mat'),'avg_w_time','-v7.3');

    % Compile final wavelets into cell array to be used in IG inversion
    %%
    n_win = length(unique(w_sum(1,:)));
    for i_vol = 1:1:i_vol_max    
        all_wavelets_freq{1,i_vol} = w_avg_freq(:,((i_vol-1)*n_win)+1:n_win*i_vol);% Separate  the frequency spectruml of different offset volumes
        all_wavelets_time{1,i_vol} = w_avg_time(:,((i_vol-1)*n_win)+1:n_win*i_vol);% Separate  the  wavelet of different offset volumes
        wavelet_length(i_vol) = size(w_avg_time,1)-2;
    end

    [~,min_ind] = min(wavelet_length);
    wavelets.min_wavelet_z_grid = all_wavelets_time{1,min_ind}(1,:);
    [~,max_ind] = max(wavelet_length);
    wavelets.max_wavelet_z_grid = all_wavelets_time{1,max_ind}(1,:);
    wavelets.all_wavelets_time = all_wavelets_time;
    wavelets.all_wavelets_freq = all_wavelets_freq;
    wavelet_save_path = strcat(job_meta.wav_directory,'all_wavelets_time.mat');
    save(wavelet_save_path,'-struct','wavelets','-v7.3');    
    save(job_meta_path,'-struct','job_meta','-v7.3');    
end
