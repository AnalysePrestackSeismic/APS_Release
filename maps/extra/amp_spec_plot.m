function [] = amp_spec_plot(pltdata,samprate)
%---------------------------------------
%
% average amp spectrum of some data
%
%--------------------------------------
%
nt = length(pltdata);
dt = samprate/1.e6;

% invariants
fsample = 1./dt;

% freq content
%figure(1)
%subplot(2,1,1)
avfft = mean(abs(fft(pltdata)),2);

figure; plot((0:round(nt/2))*fsample/nt,20*log10(avfft(1:round(nt/2)+1))  )
axis tight
xlabel('Frequency (Hz)')
ylabel('Power dB')
title('Amplitude spectrum')
