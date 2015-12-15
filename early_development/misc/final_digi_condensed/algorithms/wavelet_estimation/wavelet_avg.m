function [] = wavelet_avg(job_meta_path)
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

for i_vol = 1:1:job_meta.nvols
    for i_block = 1:1:job_meta.n_blocks 
        fid_wav = fopen(strcat(job_meta.wav_directory,job_meta.volumes{i_vol}, ...
            '_fft_wavelets_block_',num2str(i_block),'.bin'));
        w = fread(fid_wav,'float32');
        n_win = w(1);
        fclose(fid_wav);
        w = reshape(w(2:end),[],n_win);
        wav_win_ind = w(1,:);        
        if size(w) > 1
            if i_block == 1
                tmp_w = w;
            else
               
%                 if sum(logical(w(1,:))) > sum(logical(tmp_w(1,:)))                    
%                     tmp_w(1,:) = w(1,:);
%                 end
            [Loca,~] = ismember(tmp_w(1,:),wav_win_ind);
            tmp_w(2:end,Loca) = tmp_w(2:end,Loca)+w(2:end,Loca);
            end
        end
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