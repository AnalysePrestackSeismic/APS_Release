function wfilter = match(reference,shapee,prewhite,file_out)

% USE:      type: filter = match(reference,shapee,prewhite,'file_out');
%
% INPUTS:   reference = wavelet to match to
%           shapee = wavelet to shape
%           prewhite = prewhitenning as dB down from 0 power (should be a 
%                           negative number)
%           file_out = 'filename for output filter' (must be entered within
%                           quotation marks)
%
% OUTPUTS:  filter = time domain least-squares matching filter to match
%                       shapee to reference. Filter is scaled such that the
%                       RMS of shapee is the same after filtering as
%                       before.

% get number of samples of inputs
nt = length(reference);
%nt=nt*4;

% prewhiten shapee
% fa_reference = abs(fft(reference,nt));
% fp_reference = 20*log10(fa_reference/max(fa_reference));
% bandwidth = (fp_reference > prewhite);
% reference = circshift(ifft(fa_reference),floor(nt/2));
% 
% fa_shapee = abs(fft(shapee,nt));
% fp_shapee = 20*log10(fa_shapee/max(fa_shapee));
% edges = logical(abs([0;diff(bandwidth)]));
% edges_idx = cumsum(edges);
% low = (edges_idx == 0);
% mid = (edges_idx == 2);
% high = (edges_idx == 4);
% fp_shapee(low) = fp_shapee(and(high,edges));
% fp_shapee(mid) = fp_shapee(and(mid,edges));
% fp_shapee(high) = fp_shapee(and(high,edges));
% faw_shapee = (10.^(fp_shapee/20))*max(fa_shapee);
% wshapee = circshift(ifft(faw_shapee),floor(nt/2));
% 
% shapee = circshift(ifft(fa_shapee),floor(nt/2));
% 
% white = zeros(nt,1);
% white(floor((nt/2)+1),1) = max(wshapee)-max(shapee);
% wshapee = shapee+white;

wshapee = shapee;

% plot inputs
figure(1)
set(1,'Units','inches','Position',[0 0 20 15]);
subplot(2,2,1)
plot(reference,'LineWidth',3)
hold all
plot(shapee,'LineWidth',3)
axis tight

% white = zeros(nt,1);
% white(floor((nt/2)+1),1) = norm(shapee)*(prewhite/100);
% wshapee = shapee+white;

% make convolution matrix
% wG = convmtx(wshapee,nt);   % with pre-whitenning
% G = convmtx(shapee,nt);     % without pre-whitenning
clip = (((2*nt)-1)-nt)/2;
if round(clip)==clip    
%     wG = wG(clip+1:end-clip,:);
%     G = G(clip+1:end-clip,:);    
% code for when user has no signal processing toolbox (avoids convmtx)
    for i = 1:nt
        wG(i,:) = circshift(wshapee,i+clip);
        G(i,:) = circshift(shapee,i+clip);
    end
    A = spdiags(wG,(-clip:1:clip));
    wG = spdiags(A,(-clip:1:clip),nt,nt);
    B = spdiags(G,(-clip:1:clip));
    G = spdiags(B,(-clip:1:clip),nt,nt);
else 
    clip = round(clip);
%     wG = wG(clip+1:end-clip+1,:);   % with pre-whitenning
%     G = G(clip+1:end-clip+1,:);     % without pre-whitenning
% code for when user has no signal processing toolbox (avoids convmtx)    
    for i = 1:nt
        wG(i,:) = circshift(wshapee,i+clip-1);
        G(i,:) = circshift(shapee,i+clip-1);
    end
    A = spdiags(wG,(-clip:1:clip));
    wG = spdiags(A,(-clip:1:clip),nt,nt);
    B = spdiags(G,(-clip:1:clip));
    G = spdiags(B,(-clip:1:clip),nt,nt);
end

% calculate least-squares matching filter
wfilter = lsqr(wG,reference,1e-4,1000);  % with pre-whitenning
filter = lsqr(G,reference,1e-4,1000);    % without pre-whitenning

% ensure RMS of shapee remains the same after filtering
wscalar = norm(shapee)/norm(G*wfilter);     % with pre-whitenning
wfilter = wfilter*wscalar;                  % with pre-whitenning
scalar = norm(shapee)/norm(G*filter);       % without pre-whitenning
filter = filter*scalar;                     % without pre-whitenning

% plot outputs
plot(G*wfilter,'LineWidth',3,'linestyle','--')
title('Time domain shaping','FontSize',14)
ylabel('Amplitude','FontSize',14)
xlabel('Sample Number','FontSize',14)
%legend('Reference wavelet','Wavelet to be shaped','Wavelet after shaping')
legend('Reference wavelet','Wavelet to be shaped','Wavelet after shaping','Location','SouthEast')
set(gca,'FontSize',14);

subplot(2,2,3)
plot(wfilter,'LineWidth',3)
hold all
plot(filter,'red','LineWidth',3,'linestyle','--')
title('Time domain filter','FontSize',14)
ylabel('Amplitude','FontSize',14)
xlabel('Sample Number','FontSize',14)
legend('With pre-whitenning','Without pre-whitenning')
set(gca,'FontSize',14);
axis tight

% make data for frequency domain plots
f_reference = abs(fft(reference));
f_reference = f_reference(1:round(nt/2));
f_reference = 20*log10(f_reference/max(f_reference));

f_shapee = abs(fft(shapee));
f_shapee = f_shapee(1:round(nt/2));
f_shapee = 20*log10(f_shapee/max(f_shapee));

f_shaped = abs(fft(G*wfilter));
f_shaped = f_shaped(1:round(nt/2));
f_shaped = 20*log10(f_shaped/max(f_shaped));

f_wfilter = abs(fft(wfilter));
f_wfilter = f_wfilter(1:round(nt/2));
f_wfilter = 20*log10(f_wfilter/max(f_wfilter));

f_filter = abs(fft(filter));
f_filter = f_filter(1:round(nt/2));
f_filter = 20*log10(f_filter/max(f_filter));

octaves = [125,125,125/2,125/2,125/4,125/4,125/8,125/8,125/16,125/16,125/32,125/32,125/64,125/64,125/128,125/128;5,-40,5,-40,5,-40,5,-40,5,-40,5,-40,5,-40,5,-40]';

subplot(2,2,2)
plot((0:1:(round(nt/2)-1))*125/(round(nt/2)-1),f_reference,'LineWidth',3)
axis([0,125,-40,5])
hold all
plot((0:1:(round(nt/2)-1))*125/(round(nt/2)-1),f_shapee,'LineWidth',3)
plot((0:1:(round(nt/2)-1))*125/(round(nt/2)-1),f_shaped,'LineWidth',3,'linestyle','--')
for i=1:8
    plot(octaves(1+(2*(i-1)):2*i,1),octaves(1+(2*(i-1)):2*i,2),'black','linestyle','--');
end
title('Frequency domain shaping','FontSize',14)
ylabel('Power (dB below maximum amplitude)','FontSize',14)
xlabel('Frequency (Hz)','FontSize',14)
legend('Reference spectrum','Spectrum to be shaped','Spectrum after shaping')
set(gca,'FontSize',14);

subplot(2,2,4)
plot((0:1:(round(nt/2)-1))*125/(round(nt/2)-1),f_wfilter,'LineWidth',3)
axis([0,125,-40,5])
hold all
plot((0:1:(round(nt/2)-1))*125/(round(nt/2)-1),f_filter,'red','LineWidth',3,'linestyle','--')
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
print('-r300','-dtiff',sprintf('%s.tiff',file_out));

dlmwrite(sprintf('%s_with_prewhitening.txt',file_out),wfilter,'\t');
dlmwrite(sprintf('%s_without_prewhitening.txt',file_out),filter,'\t');

end
