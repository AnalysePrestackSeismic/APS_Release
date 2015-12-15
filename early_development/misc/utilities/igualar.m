function [filter] = igualar(reference,shapee,prewhite)
% USE:         Matching Wavelets via FREQUENCY Domain Division
% INPUTS:      reference = wavelet to match to
%              shapee = wavelet to shape
%              prewhite = percentage of energy of reference fft, apply as
%              fraction e.g. 0.04 as 4%
%
% IGUALAR = MATCH, in Portuguese
%

% counting samples
nt = length(reference)
nts = length(shapee)

% make array of zeros for padding
pad = zeros(floor((abs(length(reference)-length(shapee)))/2),1);

% pad the shorter wavelet
if length(reference) > length(shapee)
          shapee = [pad;shapee;pad];
else
          reference = [pad;reference;pad];
end

% padded ref - for freq domain plots
ntpr = length(reference)

% Prewhitening - via % of energy of reference wavelet (fft)
pw = prewhite
%e = sum(abs(reference).^2) 
e = sum((abs(fft(reference))).^2)./nt
pw = pw*e

% make the matching operator by frequency domain division to match wavelet 2 to wavelet 1
% matches amplitude and phase spectra
filter_fft = fft(conv(reference,shapee))./(fft(conv(shapee,shapee))+pw);
filter_fft_nopw = fft(conv(reference,shapee))./fft(conv(shapee,shapee));

% make the matching operator by frequency domain division to match wavelet 2 to wavelet 1
% only matches amplitude spectra
% filter_fft = abs(fft(conv(reference,shapee)))./(abs(fft(conv(shapee,shapee)))+pw);

% Estimating shift
if nts > nt
          shift = floor(((2*nts)-1)/2)
else 
          shift = floor(((2*nt)-1)/2)
end

% apply the operator by time domain convolution to make wavelet 3 (wavelet 3 ~= wavelet 1)
shaped = conv(shapee,circshift(ifft(filter_fft),shift),'same');

filter = circshift(ifft(filter_fft),shift)
filter_nopw = circshift(ifft(filter_fft_nopw),shift)

% plot outputs
figure(1)
set(1,'Units','inches','Position',[0 0 20 15]);
subplot(2,2,1)
plot(reference,'LineWidth',3)
hold all
plot(shapee,'LineWidth',3);
axis tight

plot(shaped,'LineWidth',3,'linestyle','--')
axis tight
title('Time domain shaping','FontSize',14)
ylabel('Amplitude','FontSize',14)
xlabel('Sample Number','FontSize',14)
legend('Reference wavelet','Wavelet to be shaped','Wavelet after shaping')
set(gca,'FontSize',14);

subplot(2,2,3)
plot(filter,'LineWidth',3);
hold all
plot(filter_nopw,'cyan','LineWidth',3,'linestyle','--')
title('Time domain filter','FontSize',14)
ylabel('Amplitude','FontSize',14)
xlabel('Sample Number','FontSize',14)
legend('With pre-whitenning','Without pre-whitenning')
set(gca,'FontSize',14);
axis tight

% make data for frequency domain plots
f_reference = abs(fft(reference));
f_reference = f_reference(1:round(ntpr/2));
f_reference = 20*log10(f_reference/max(f_reference));

f_shapee = abs(fft(shapee));
f_shapee = f_shapee(1:round(ntpr/2));
f_shapee = 20*log10(f_shapee/max(f_shapee));

f_shaped = abs(fft(shaped));
f_shaped = f_shaped(1:round(ntpr/2));
f_shaped = 20*log10(f_shaped/max(f_shaped));

f_filter = abs(fft(filter));
f_filter = f_filter(1:round(ntpr/2));
f_filter = 20*log10(f_filter/max(f_filter));

f_filter_nopw = abs(fft(filter_nopw));
f_filter_nopw = f_filter_nopw(1:round(ntpr/2));
f_filter_nopw = 20*log10(f_filter_nopw/max(f_filter_nopw));

octaves = [125,125,125/2,125/2,125/4,125/4,125/8,125/8,125/16,125/16,125/32,125/32,125/64,125/64,125/128,125/128;5,-40,5,-40,5,-40,5,-40,5,-40,5,-40,5,-40,5,-40]';

subplot(2,2,2)
plot((0:1:(round(ntpr/2)-1))*125/(round(ntpr/2)-1),f_reference,'LineWidth',3)
axis([0,125,-40,5])
hold all
plot((0:1:(round(ntpr/2)-1))*125/(round(ntpr/2)-1),f_shapee,'LineWidth',3)
plot((0:1:(round(ntpr/2)-1))*125/(round(ntpr/2)-1),f_shaped,'LineWidth',3,'linestyle','--')
for i=1:8
    plot(octaves(1+(2*(i-1)):2*i,1),octaves(1+(2*(i-1)):2*i,2),'black','linestyle','--');
end
title('Frequency domain shaping','FontSize',14)
ylabel('Power (dB below maximum amplitude)','FontSize',14)
xlabel('Frequency (Hz)','FontSize',14)
legend('Reference spectrum','Spectrum to be shaped','Spectrum after shaping')
set(gca,'FontSize',14);

subplot(2,2,4)
plot((0:1:(round(ntpr/2)-1))*125/(round(ntpr/2)-1),f_filter,'LineWidth',3)
axis([0,125,-40,5])
hold all
plot((0:1:(round(ntpr/2)-1))*125/(round(ntpr/2)-1),f_filter_nopw,'cyan','LineWidth',3,'linestyle','--')
for i=1:8
    plot(octaves(1+(2*(i-1)):2*i,1),octaves(1+(2*(i-1)):2*i,2),'black','linestyle','--');
end
title('Frequency domain filter','FontSize',14)
ylabel('Power (dB below maximum amplitude)','FontSize',14)
xlabel('Frequency (Hz)','FontSize',14)
legend('With pre-whitenning','Without pre-whitenning','Location','SouthEast')
set(gca,'FontSize',14);
set(gcf, 'PaperUnits', 'inches');                               %set dimensions to inches
set(gcf, 'PaperPosition', [1,1,12,9]);                         %set picture size
% print('-r300','-dtiff',sprintf('%s.tiff',file_out));

% dlmwrite(sprintf('%s_with_prewhitening.txt',file_out),wfilter,'\t');
end