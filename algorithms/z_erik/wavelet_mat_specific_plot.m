function wavelet_mat_specific_plot(job_meta_path,wave_mat_path)
%Function to plot wavelets in the formatted .mat files
%
%%
% Load job meta information 
job_meta = load(job_meta_path);
wave = load(wave_mat_path);    
n_vol = size(wave.all_wavelets_freq,2);
[ns,n_win]=size((wave.all_wavelets_freq{1,2}));
n_s=ns-1;
w_freq=zeros(n_s,n_win,n_vol);
w_time=zeros(n_s,n_win,n_vol);

start_vol=1;
end_vol=6;

for i_vol=1:end_vol
    w_freq(:,:,i_vol) = wave.all_wavelets_freq{1,i_vol}(2:end,:);

    w_time(:,:,i_vol) = wave.all_wavelets_time{1,i_vol}(2:end,:);
end

plotinc = 1;
maxplot = end_vol;
startplot = start_vol;
totalplots = ((maxplot - startplot)/plotinc ) + 1;

freq_axis = (1e6/job_meta.s_rate)/2*linspace(0,1,job_meta.ns_win/2);
figure(91)
for i_vol = startplot:plotinc:maxplot
%title(i_block)
subplot(3,totalplots,(((i_vol-startplot)/plotinc)+1)); imagesc(w_time(:,:,i_vol)); xlim([1 14.5]);
subplot(3,totalplots,(((i_vol-startplot)/plotinc)+1)+totalplots); plot(freq_axis,w_freq(3:2+ceil(job_meta.ns_win/2),:,i_vol)); set(gca,'XTick',1:ceil(10*freq_axis(2)):ceil(freq_axis(end))); xlim([1 100])
subplot(3,totalplots,((((i_vol-startplot)/plotinc)+1)+(totalplots*2))); imagesc(w_freq(3:2+ceil(job_meta.ns_win/2),:,i_vol)'); ylim([0.5 16.5])
end

figure(92)
for i_vol = startplot:plotinc:maxplot
    subplot(6,ceil(totalplots/6),((i_vol-startplot)/plotinc)+1); plot(w_time(:,:,i_vol)); set(gca,'ActivePositionProperty', 'Position')
end
end

