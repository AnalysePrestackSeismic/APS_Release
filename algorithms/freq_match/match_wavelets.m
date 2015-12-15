function match_wavelets


% constants ===============================================================
NF = 91; %matching filter length
percenttaper = 5;  % taper to apply to the start and end of wavelet
                   %e.g 5 is 5% of the length of the wavelet at the start and at the end

%==========================================================================

load('wavelets_for_james.mat');

wref = wavelets.all_wavelets_time{1,9}(:,2);
wout = wavelets.all_wavelets_time{1,33}(:,2);



if floor(NF/2) == NF/2;
    NF = NF + 1;
end





wref = taperwavelet(wref,percenttaper);
wout = taperwavelet(wout,percenttaper);

[filter,output] = calculate_filter(wref,wout,NF,0);

%avfftfilt = amp_spec(filter,4);
%avfftwref = amp_spec(wref,4);
%avfftoutput = amp_spec(output,4);

figure(2)
axis tight
subplot(3,2,1); plot(wref)
title('input wavelet')
subplot(3,2,2); plot(wout)
title('desired wavelet')
subplot(3,2,3); plot(filter)
title('weiner filter')
subplot(3,2,4); plot(wref,'-b');
title('input wavelet in blue, desired in green, matched in red')
hold all; plot(wout,'-g'); plot(output,'--r');
hold off
subplot(3,2,5); amp_spec(filter,4);
xlabel('Frequency (Hz)')
ylabel('Power dB')
title('Weiner filter amp spectrum')
subplot(3,2,6); amp_spec(wref,4); 
hold all; 
amp_spec(output,4);
hold off
xlabel('Frequency (Hz)')
ylabel('Power dB')
title('Amplitude spectrum, input blue, output green')



end

function [] = amp_spec(pltdata,samprate)
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
pwrspec = mean(abs(fft(pltdata)),2);
plot((0:round(nt/2))*fsample/nt,20*log10(pwrspec(1:round(nt/2)+1))  )
end

function [taper_wav] = taperwavelet(wavelet,percenttaper)
% make taper to apply to signal before fft
taperlen = floor((length(wavelet)*0.01)*percenttaper);
%taperst = linspace(0,1,taperlen)';
taperst = (sin(linspace((-pi/2),(pi/2),taperlen)')+1)/2;
taperend = 1 - taperst;
taperapply = [taperst;ones((length(wavelet)-(taperlen*2)),1);taperend];
taper_wav = bsxfun(@times,wavelet,taperapply);
end

function [filter,output] = calculate_filter(wref,wout,NF,mu) 
  % wref wout
%mu:  prewhitening as percentage of the zero lag autocorrelation

 % Compute auto-correlation of reference wavelet
 % this was an auto correlation using toolbox
 %autoreft = xcorr(wref,wref,floor(NF/2)); 

 % this is doing it in the frequency domain
 corrLength=length(wref)+length(wref)-1;
 midpos = length(wref);
 
 autorefall = fftshift(ifft(fft(wref,corrLength).*conj(fft(wref,corrLength)))); 
 autoref = autorefall((midpos - 2 - (floor(NF/2)-2)):(midpos + (floor(NF/2))));

 
 autoref = spdiags(repmat(autoref,1,NF),[-floor(NF/2):1:floor(NF/2)]); 
  
 % add pre-whitening
 autoref = autoref+((autoref(1,1)*mu/100)*eye(NF)); 
 
 % Compute cross-correlation of reference wavelet with desired wavelet
 %crossref = xcorr(wref,wout,floor(NF/2));
 
 % Compute cross-correlation of reference wavelet with desired wavelet
 corrLength=length(wref)+length(wout)-1;
 crossref = fftshift(ifft(fft(wref,corrLength).*conj(fft(wout,corrLength))));  
 %crossrefsub = crossref((midpos - 2 - (floor(NF/2)-2)):(midpos + (floor(NF/2))));
 
 %tol = 1e-8;
 %axit = 100;
 %filter = lsqr(autoref,crossref,tol,maxit);   
 filter = autoref\crossref((midpos - 2 - (floor(NF/2)-2)):(midpos + (floor(NF/2))));
 
 output = conv(wref,filter,'same'); 
 
 
end


function [filter,output] = calculate_filter_withplots(wref,wout,NF,mu) 
  % wref wout
%mu:  prewhitening as percentage of the zero lag autocorrelation

 NW = length(wref);              % length of the input
 
 if floor(NF/2) == NF/2;
     NF = NF + 1;
 end
  
 % Compute auto-correlation of reference wavelet
 % Cref = convmtx(wref,NF);
 % autoref3 =  Cref'*Cref; 
 
 figure(9)
 
 % this was an auto correlation using toolbox
 %autoreft = xcorr(wref,wref,floor(NF/2)); 
 subplot(3,1,1); plot(autoreft); title('autoreft')
 % this is doing it in the frequency domain
 corrLength=length(wref)+length(wref)-1;
 midpos = length(wref);
 autorefall = fftshift(ifft(fft(wref,corrLength).*conj(fft(wref,corrLength)))); 
 autoref = autorefall((midpos - 2 - (floor(NF/2)-2)):(midpos + (floor(NF/2))));
 autoreforig = autoref;
 subplot(3,1,2); plot(autoref); title('autoref')
 subplot(3,1,3); plot(autoreft); title('autoref')
 hold on
 %subplot(3,1,3); plot([zeros((midpos - 2 - (floor(NF/2)-1)),1); autoreft],'-r'); title('autoreft')
 subplot(3,1,3); plot(autoref,'r'); title('autoref')
 hold off
 figure(1)
 subplot(2,3,1); plot(autoref); title('autoref')
 
 autoref = spdiags(repmat(autoref,1,NF),[-floor(NF/2):1:floor(NF/2)]); 
 subplot(2,3,2); imagesc(autoref); title('autoref operator')
 subplot(2,3,3); plot(autoref(1,:)); title('row of autoref operator')
  
 % add pre-whitening
 autoref = autoref+((autoref(1,1)*mu/100)*eye(NF)); 
 subplot(2,3,4); imagesc(autoref); title('autoref operator with pre-whitening')
 subplot(2,3,5); plot(autoref(1,:)); title('row of autoref operator with pre-whitening')
 
 % Compute cross-correlation of reference wavelet with desired wavelet
 %crossref = xcorr(wref,wout,floor(NF/2));
 % Compute cross-correlation of reference wavelet with desired wavelet
 corrLength=length(wref)+length(wout)-1;
 crossref = fftshift(ifft(fft(wref,corrLength).*conj(fft(wout,corrLength))));  
 crossrefsub = crossref((midpos - 2 - (floor(NF/2)-2)):(midpos + (floor(NF/2))));
 
 %pad = NW+NF-1-length(wout);
 %crossref = Cref'*[wout; zeros(pad,1)];
 
 subplot(2,3,6); plot(autoreforig); title('autoreforig')
 hold on
 subplot(2,3,6); plot(crossrefsub,'-r'); title('crossrefsub')
 hold off
 %tol = 1e-8;
 %axit = 100;
 %filter = lsqr(autoref,crossref,tol,maxit);   
 filter = autoref\crossrefsub;
 
 output = conv(wref,filter,'same'); 
 
 
end
                              



