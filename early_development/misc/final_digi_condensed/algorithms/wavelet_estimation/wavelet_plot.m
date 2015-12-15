function wavelet_plot(job_meta_path,i_block)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Load job meta information 
job_meta = load(job_meta_path);

for i_vol = 1:1:job_meta.nvols
    fid_wav = fopen(strcat(job_meta.wav_directory,job_meta.volumes{i_vol},'_fft_wavelets_block_',i_block,'.bin'));
    w = fread(fid_wav,'float32');
    w_freq{i_vol} = reshape(w,[],job_meta.n_win);
    % w_freq{i_vol} = w(3:end,:);
    % Convert from frequency domain to time domain
    w_time{i_vol} = circshift(ifft(w_freq{i_vol}(3:end,:),'symmetric'),floor(job_meta.ns_win/2));
end

%imagesc(w_time{1});

figure
subplot(1,3,1); imagesc(cell2mat(w_time')');
%colorbar
%colormap('gray')

freq_axis = (1e6/job_meta.s_rate)/2*linspace(0,1,job_meta.ns_win/2+1);

subplot(1,3,2); plot(freq_axis,w(3:2+ceil(job_meta.ns_win/2),:))

subplot(1,3,3); imagesc(w_freq{1}(3:2+ceil(job_meta.ns_win/2),:)')

end

