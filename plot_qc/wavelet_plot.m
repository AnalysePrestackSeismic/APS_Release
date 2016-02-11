function wavelet_plot(job_meta_path,i_block)
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
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


% Load job meta information 
job_meta = load(job_meta_path);


%for i_vol = 1:1:job_meta.nvols
    fid_wav = fopen(strcat(job_meta.wav_directory,'fft_wavelets_block_',i_block,'.bin'));
    w_tmp = fread(fid_wav,'float32');
    fclose(fid_wav);
    %n_win = w_tmp(1);
    n_vol = w_tmp(1);
    n_win = w_tmp(2);
    stdev = w_tmp(3); 
    live_offset = w_tmp(4);
    wb_z_avg = w_tmp(5);
    %w_freq{i_vol} = reshape(w,[],n_win);
    %w_freq{i_vol} = reshape(w_tmp(3:end),[],n_vol);
    wave_offset = 6;
    w_freq = reshape(w_tmp(wave_offset:end),[],n_win,n_vol);
    % w_freq{i_vol} = w(3:end,:);
    % Convert from frequency domain to time domain
    w_time = circshift(ifft(w_freq(wave_offset:end,:,:),'symmetric'),floor(job_meta.ns_win/2));
    %w_time = circshift(ifft(w_freq(3:end,:,:),'symmetric'),floor(size( ifft(w_freq(3:end,:,:),'symmetric'),1)/2));
    %w_time = ifft(w_freq(3:end,:,:),'symmetric');
%end


plotinc  = 1;
maxplot = n_vol;
%maxplot = 15;
startplot = 1;
totalplots = ((maxplot - startplot)/plotinc ) + 1;



%imagesc(w_time{1});
freq_axis = (1e6/job_meta.s_rate)/2*linspace(0,1,job_meta.ns_win/2);
figure(91)
for i_vol = startplot:plotinc:maxplot
title(i_vol) % i_block
subplot(3,totalplots,(((i_vol-startplot)/plotinc)+1)); imagesc(w_time(:,:,i_vol));
subplot(3,totalplots,(((i_vol-startplot)/plotinc)+1)+totalplots); plot(freq_axis,w_freq(3:2+ceil(job_meta.ns_win/2),:,i_vol)); set(gca,'XTick',1:ceil(10*freq_axis(2)):ceil(freq_axis(end)))
subplot(3,totalplots,((((i_vol-startplot)/plotinc)+1)+(totalplots*2))); imagesc(w_freq(3:2+ceil(job_meta.ns_win/2),:,i_vol)');
end
%colorbar
%colormap('gray')
figure(92)
for i_vol = startplot:plotinc:maxplot
    subplot(6,ceil(totalplots/6),((i_vol-startplot)/plotinc)+1); plot(w_time(:,:,i_vol)); set(gca,'ActivePositionProperty', 'Position')
end

% xmin = 1;
% xmax = 129;
% wavelets_vec = cell2mat(wavelet_norm);
% wavelets_vec = wavelets_vec(:);
% ymin = min(wavelets_vec);
% ymax = max(wavelets_vec);


% for ii = 1:1:ns
%     
%    subplot(1,totalvol,1); plot(wavelet_norm{1}(:,ii)); axis([xmin xmax ymin ymax])
%    subplot(1,totalvol,2); plot(wavelet_norm{2}(:,ii)); axis([xmin xmax ymin ymax])
%    subplot(1,totalvol,3); plot(wavelet_norm{3}(:,ii)); axis([xmin xmax ymin ymax])
%    subplot(1,totalvol,4); plot(wavelet_norm{4}(:,ii)); axis([xmin xmax ymin ymax])
%    subplot(1,totalvol,5); plot(wavelet_norm{5}(:,ii)); axis([xmin xmax ymin ymax])
%    title(num2str(ii))
%    pause(0.1)
% end





end

