function wavelet_mat_specific_plot(job_meta_path,wave_mat_path)
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

