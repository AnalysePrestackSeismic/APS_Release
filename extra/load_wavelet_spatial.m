function [ wavelets ] = load_wavelet_spatial(job_meta_path,wavelet_dir_path,i_block,algo)
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
%   Function loads Spatially Varying wavelet for ith block using a chosen
%   algorithm by user
% Input:
% algo passed as a string
% i_block : Block Number Passed as Number

% Output:
% Algo: 1: For using Raw wavelets
% Algo 2:  FOr arithmetic Average of Wavelets in Defined Neighbourhood
%----------------------------------------
%%
job_meta = load(job_meta_path);             % Load job meta information
if job_meta.is_gather == 0                  % For the case of using angle stacks
    i_vol_max = job_meta.nvols;             % Number of angle stacks
else                                        % For the case of using angle gathers
    i_vol_max = ((job_meta.tkey_max -  job_meta.tkey_min )/ job_meta.tkey_inc) +  1; % Number of angles in angle gather
end
loopfin = size(job_meta.liveblocks,1);          % Number of live blocks
%%
%Creating Look up Table
ils = sort(unique(job_meta.block_keys(:,1)));   % Sort inlines of  fist corner point of blocks in increasing order and remove duplicates
xls = sort(unique(job_meta.block_keys(:,3)));   % Sort xlines of  fist corner point of blocks in increasing order and remove duplicates

ilslen = length(ils);                           % Number of unique inlines of  fist corner point of blocks
xlslen = length(xls);                           % Number of unique inlines of  fist corner point of blocks
ilxl_lookup = zeros(ilslen,xlslen,2);           % Create look up table for inlines and cross-lines
lpi = 1;                                        % Loop Index
%Loop for creating look up table from Job Meta Information
while lpi <= loopfin
    cur_block = job_meta.liveblocks(lpi);                                 % Block Number for Current Live Block
    currow = find(ils == job_meta.block_keys(cur_block,1),1,'first');     % Find ROW = The inline matching ith block first coner point inline
    curcol = find(xls == job_meta.block_keys(cur_block,3),1,'first');     % Find COULMN = The inline matching ith block first coner point xline
    
    ilxl_lookup(currow,curcol,1) = cur_block;                             % Store Block Number in the first cell indexed by found ROW and COLUMN
    ilxl_lookup(currow,curcol,2) = job_meta.stdev(cur_block);             % Store Block wavelet variance in the 2nd cell indexed by found ROW and COLUMN
    ilxl_lookup(currow,curcol,3) = job_meta.wb_z_avg(cur_block);          % Store Block Average Water Bottom in the 3rd cell indexed by found ROW and COLUMN
    ilxl_lookup(currow,curcol,4) = job_meta.live_offset(cur_block);       % Store Block Live Offset Value  in the 4th cell indexed by found ROW and COLUMN
    lpi = lpi + 1;                                                        % Increment Loop Index
end
Row_max = size(ilxl_lookup,1);                                            % Max number of rows in Block Look up matrix
Col_max = size(ilxl_lookup,2);                                            % Max number of columns in Block Look up matrix
%--------Plotting-------------
% figure(1); imagesc(ilxl_lookup(:,:,2)); title('input variance');%caxis([500 1000]);                       % Plot Variance
% figure(11); imagesc(ilxl_lookup(:,:,3)); title('input water bottom depth'); %caxis([500 1000]);           % Plot Average Water Bottom Z
% figure(12); imagesc(ilxl_lookup(:,:,4)); title('input live offset'); %caxis([0 50]);                      % Plot Live Offset
% figure(14); imagesc(ilxl_lookup(:,:,1)); title('block number'); %caxis([0 50]);                           % Plot Block Number
clear ils xls ilslen xlslen lpi cur_block currow curcol;                                                    %Clear redundant Variables

%%
% Switch Case fo choosing algorithm of wavelet smoothening
% Case 1 : Raw discreet wavelet for individual Block used
% Case 2 : Smoothened wavelet arithmatic averaged over the adjacent live neighbours
% Case 3 : ...
switch algo
    case '1'
        s= strcat(wavelet_dir_path,'fft_wavelets_block_',num2str(i_block),'.bin');
        fid_wav = fopen(s);
        w = fread(fid_wav,'float32');                                       % Load the wavelet for ith block
        n_vol = w(1);               % Number of volumes (angle stacks)
        n_win = w(2);               % Number of windows for wavelet estimation
        stdev = w(3);               % Standard Dev:  might like to think about calculating this properly based on live samples
        live_offset = w(4);         % The live angle volume
        wb_z_avg = w(5);            % Average Z waterbottom in the block?
        wave_offset = 6;            % This is the index in the wavlet file from where the individual wavelet data starts
        neig_rad = 0;                                                % Define neighbourhood Radius
        w_proc = zeros(n_win*(job_meta.ns_win+2),n_vol);             % Matrix for storing all wavelets for all volumes (Note: length of wavlet array = number of samples in wavelet) + 2
        w_proc = reshape(w(wave_offset:end),[],n_vol);               % Reshape the Matrix sperating the different angle volumes in different columns
        w_all = w_proc;
    
    case '2'
        i_block=str2double(i_block);
        [i_block_row,i_block_col] = find (ilxl_lookup(:,:,1)==i_block);         % Find row and column of ith block in lookup table
        neig_rad = 2;                                                           % Define neighbourhood Radius
        no_neig = (1+2*neig_rad)^2;
        neig_blocks = zeros(no_neig,1);                                         % Matrix for storing neighbourhood blocks
                
        %Loop to find the block numbers of neighbouring blocks
        count_b = 0;                                                            % Initialize Counter for neighbouring blocks
        
        for rb= (i_block_row-neig_rad):(i_block_row+neig_rad)                   % Loop through rows
            for cb= (i_block_col-neig_rad):(i_block_col+neig_rad)               % Loop though columns
                if(rb>0 && cb>0 && rb <=Row_max &&cb <=Col_max)                      % Check if the hypothetical neghbouring block is outsize the outer boundary
                    count_b = count_b+1;                                        % Count index (counter) of valid live neighbouring block
                    neig_blocks(count_b)=ilxl_lookup(rb,cb,1);                  % Store block number of neighbour
                end
            end
        end
        neig_blocks = neig_blocks(neig_blocks~=0);                              % Leave non-live Blocks out
        no_neig = length(neig_blocks);                                           %Total number of live neighbours
        clear rb cb count_b;
        
        % Loop tthrough Neighbouring Blocks for smoothening
        count = 1;
        for lpi=1: length(neig_blocks)
            s= strcat(wavelet_dir_path,'fft_wavelets_block_',num2str(neig_blocks(lpi)),'.bin');
            s_exist = exist(s,'file');                                  % check to see that the wavelet file actualy exists
            if s_exist > 0
                fid_wav = fopen(s);                                     % Open the binary wavelet file for this block
                w = fread(fid_wav,'float32');                           % Read the binary file and convert to 32 bit floating point precission
                n_vol = w(1);                                           % Number of volumes (angle stacks)
                n_win = w(2);                                           % Number of windows for wavelet estimation
                %             stdev = w(3);                                           % Standard Deviation :might like to think about calculating this properly based on live samples
                %             live_offset = w(4);                                     % The live angle volume
                %             wb_z_avg = w(5);                                        % Average Z waterbottom in the block?
                wave_offset = 6;                                        % This is the index in the wavlet file from where the individual wavelet data starts
                fclose(fid_wav);                                        % Close file
                if count == 1                                           % Initialize the wavelet matrix if this is the first live block
                    w_all = zeros(n_win*(job_meta.ns_win+2),n_vol,length(neig_blocks)); % Matrix for storing all wavelets for all volumes for all live blocks, (Note: length of wavlet array = number of samples in wavelet) + 2
                    w_all(:,:,count) = reshape(w(wave_offset:end),[],n_vol);% Reshape the Matrix sperating the different angle volumes in different columns for the first block
                    max_n_win = n_win;                                  % Intitialize maximum number of wavelet windows in block as number of windows in first block
                    max_n_vol = n_vol;                                  % Intitialize maximum number of angle volumes in block as number of windows in first block
                else                                                    % If this is not the first live block, append already initialized wavelet matrix
                    if n_win > max_n_win                                % If you encounter a block with bigger window than what has been encountered ljust append a slab for the extra window
                        w_append = zeros((n_win-max_n_win)*(job_meta.ns_win+2),n_vol,no_neig); % Create new slab of the extra length of window for all volumes for all live blocks
                        w_all = [w_all; w_append];                      % Append the new slab for sccomodating the extra length of window
                        w_all(1:n_win*(job_meta.ns_win+2),:,count) = reshape(w(wave_offset:end),[],n_vol); % Reshape the Matrix sperating the different angle volumes in different columns for the current block
                        max_n_win = n_win;
                    else                                                % If this block has the same or less number of wavelet windows than encountered before
                        w_all(1:n_win*(job_meta.ns_win+2),:,count) = reshape(w(wave_offset:end),[],n_vol);% Reshape the Matrix sperating the different angle volumes in different columns for the current block
                    end
                end
                count  = count + 1;                                     % Increment counter
            end
        end
        clear s lpi count w;
        w_proc = sum(w_all,3);                                       % Sum the wavelet spectrum including number of live traces in the blocks for all neighboring blocks
        
    otherwise
        fprintf ('Error:')
end
%%
                 
        get_windows = squeeze(w_all(1:job_meta.ns_win+2:end,1,:));  % Create a matrix of window starts
        get_windows = unique(get_windows(get_windows ~= 0));         % Create an array of unique window starts and leave out extra zeros if any
        
        w_proc(1:job_meta.ns_win+2:size(w_proc,1),:) = repmat(get_windows,1,i_vol_max);
        w_proc = reshape(w_proc,job_meta.ns_win+2,[]); % Making a matrix of window defination ( ns_win +2) vs windows
        
        
        w_avg = [w_proc(2,:); bsxfun(@rdivide,w_proc(3:end,:),w_proc(2,:))];                                % Averaging over all traces: Divide the summation of spectrums by total number of live traces  for all windows
        w_avg_time = [w_proc(1,:); circshift(ifft(w_avg(2:end,:),'symmetric'),floor(job_meta.ns_win/2))];   % Inverse fourier transform to get the time waveletl
        w_avg_freq = [w_proc(1,:); w_avg(2:end,:)];                                                         % Matrix with the the start of wavelet windows and average frequency spectrum per window
        
        n_win = length(unique(w_proc(1,:)));                                                                % Number of wavelet windows
        
        for i_vol = 1:1:i_vol_max
            wavelet_freq{1,i_vol} = w_avg_freq(:,((i_vol-1)*n_win)+1:n_win*i_vol);                          % Separate  the frequency spectruml of different offset volumes
            wavelet_time{1,i_vol} = w_avg_time(:,((i_vol-1)*n_win)+1:n_win*i_vol);                          % Separate  the  wavelet of different offset volumes
            wavelet_length(i_vol) = size(w_avg_time,1)-2;
        end
        
        % Making the wavelet Structure
        [~,min_ind] = min(wavelet_length);
        wavelets.min_wavelet_z_grid = wavelet_time{1,min_ind}(1,:);
        [~,max_ind] = max(wavelet_length);
        wavelets.max_wavelet_z_grid = wavelet_time{1,max_ind}(1,:);
        wavelets.all_wavelets_time = wavelet_time;
        wavelets.all_wavelets_freq = wavelet_freq;
        
        
        wavelet_save_path = strcat(job_meta.wav_directory,'smo_wavelet_','algo_',algo,'_',num2str(neig_rad),'_block_',num2str(i_block),'.mat');
        save(wavelet_save_path,'-struct','wavelets','-v7.3');    
        %save(job_meta_path,'-struct','job_meta','-v7.3'); 
        
end


