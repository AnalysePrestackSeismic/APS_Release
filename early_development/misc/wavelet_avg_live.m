function [] = wavelet_avg_live(job_meta_path)
% -------------------------------------------------------------------------
% WAVELET_AVG: function to create wavelet set for use in DIGI. Currently
% wavelets vary with time and angle. Spatially variant wavelets not yet
% implimented.
%   Inputs:
%       seismic_mat_path = path to .mat file created using
%       segy_make_structure function.
%       n_blocks = number of blocks to divide processing into.
%   Outputs:
%       all_wavelets_time.mat = mat file to be input into DIGI.
% -------------------------------------------------------------------------

% Load job meta information 
job_meta = load(job_meta_path);

loopfin = size(job_meta.liveblocks,1);
lpi = 1;
%array to hold the index of time windows
maxwinloclookup = 200;
winloclookup = zeros(maxwinloclookup,1);



for i_vol = 1:1:job_meta.nvols
    %wavelet_cell{i_vol} = zeros(job_meta.ns_win+2,
    
    while lpi <= 20
    %while lpi <= loopfin
    i_block = job_meta.liveblocks(lpi);
    %for i_block = 1:1:job_meta.n_blocks
        fid_wav = fopen(strcat(job_meta.wav_directory,job_meta.volumes{i_vol}, ...
            '_fft_wavelets_block_',num2str(i_block),'.bin'));
        strcat(job_meta.wav_directory,job_meta.volumes{i_vol}, ...
            '_fft_wavelets_block_',num2str(i_block),'.bin')
        
        
        w = fread(fid_wav,'float32');
        % test if n_win is the same as before
        n_vol = w(1);
        n_win = w(2);
       
       
        fclose(fid_wav);
     
        
        
        if n_vol == 0
            
        else
            % reshape to individual volumes and then second reshape to
            % individual windows
            %w = reshape(w(3:end),[],n_vol);
            %w = reshape(w(3:end),[],n_win,n_vol);
            % this assumes sequential window numbering
            if lpi == 1
                w_all = zeros(n_win*(job_meta.ns_win+2),n_vol,loopfin);
                w_all(:,:,lpi) = reshape(w(3:end),[],n_vol);              
                max_n_win = n_win;
                max_n_vol = n_vol;
            else
                if n_win > max_n_win
                    % need to append
                    w_append = zeros((n_win-max_n_win)*(job_meta.ns_win+2),n_vol,loopfin);
                    w_all = [w_all; w_append];
                    w_all(1:n_win*(job_meta.ns_win+2),:,lpi) = reshape(w(3:end),[],n_vol);
                    max_n_win = n_win;
                else
                    % write into array
                    w_all(1:n_win*(job_meta.ns_win+2),:,lpi) = reshape(w(3:end),[],n_vol);
                end
            end
                    
                %wavelet_cell{1:n_win,1:n_vol}(:,lpi) = squeeze(mat2cell(w,130,ones(n_win,1),ones(n_vol,1)));
%             else      
%                 win_logic = winloclookup(1:n_win,1) == w(1,:,1)';
%                 if sum(win_logic) == n_win;
%                     % then carry on coding
%                     wavelet_cell{1:n_win,1:n_vol}
%                     wavelet_cell{10,4}(ns_win+2,i_block) = w(:,10,4);
%                 else
%                     winloclookup(win_logic==0,1) = w(1,win_logic==0,1); 
%                     % need to append
%                 end
%             end
            %load to maxtirx for each vol and for each time window
%             for i_win = 1:1:n_win
%                 cur_win = w(1,i_win,1);
%                 %test to see if the window value exists in the winloclookup
%                 % if not make a new cell array for it
%                 for iwl = 1:1:
%                
%             end
            %n_vol = w(1,:);
            %n_win = w(2,:);
            %w = reshape(w(3:end,:),[],nvols,n_win);
            wav_win_ind = w(1,:);
            if size(w) > 1
                if lpi == 1
                    tmp_w = w;
                else
                    
                    %                 if sum(logical(w(1,:))) > sum(logical(tmp_w(1,:)))
                    %                     tmp_w(1,:) = w(1,:);
                    %                 end
                    %check to make sure that windows are the same and only
                    %take those which are the same
                    [Loca,~] = ismember(tmp_w(1,:),wav_win_ind);
                    
                    tmp_w(2:end,Loca) = tmp_w(2:end,Loca)+w(2:end,Loca);
                end
            end
        end
    lpi = lpi + 1;    
    end
        
    % Average wavelets across blocks
    avg_w = [tmp_w(1,:); bsxfun(@rdivide,tmp_w(3:end,:),tmp_w(2,:))];
    avg_w = avg_w(:,logical(1-logical(sum(isnan(avg_w)))));

    avg_w = avg_w(:,logical(sum(avg_w(2:end,:))));
    avg_w_time = [avg_w(1,:); circshift(ifft(avg_w(2:end,:),'symmetric'),floor(job_meta.ns_win/2))];

    % Save average wavelets
    save(strcat(job_meta.wav_directory,job_meta.volumes{i_vol},'_avg_w_freq.mat'),'avg_w','-v7.3');
    save(strcat(job_meta.wav_directory,job_meta.volumes{i_vol},'_avg_w_freq.mat'),'avg_w_time','-v7.3');

    % Compile final wavelets into cell array to be used in IG inversion
    all_wavelets_freq{1,i_vol} = avg_w;
    all_wavelets_time{1,i_vol} = avg_w_time;

    % Save final wavelet files
    wavelet_length(i_vol) = size(avg_w_time,2);
end

    [~,min_ind] = min(wavelet_length);
    wavelets.min_wavelet_z_grid = all_wavelets_time{1,min_ind}(1,:);
    [~,max_ind] = max(wavelet_length);
    wavelets.max_wavelet_z_grid = all_wavelets_time{1,max_ind}(1,:);
    wavelets.all_wavelets_time = all_wavelets_time;
    wavelets.all_wavelets_freq = all_wavelets_freq;
    wavelet_save_path = strcat(job_meta.wav_directory,'all_wavelets_time.mat');
    save(wavelet_save_path,'-struct','wavelets','-v7.3');    

end