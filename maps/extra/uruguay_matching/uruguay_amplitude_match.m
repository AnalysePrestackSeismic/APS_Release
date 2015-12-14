function [] = uruguay_amplitude_match
%Design frequency-matching filters
%
%Measure freqency content of an input stack and design a zero-phase filter
%to match to the specified output spectrum
%
% stack_meta_path:          meta data for input stack
% i_block:                  block number to analyse (0 = analyse all blocks)
% output path:              path for output segy dataset
% start_win:                start window (ms)
% end_win:                  end window (ms)
% freq_pairs:               freq, amplitude pairs for the output spectrum
% max_lag:                  max autocorrelation lead/lag

% Workflow as follows
% 2. Read in average autocorrelations
% 3. Taper the average autocorrelation with a cosine taper
% 4. Take square root of amplitude spectrum
% 6. Divide output spectrum by input spectrum to get filter spectrum
% 7. Inverse fft to get filter response

[input_str ilxl input_acors] = segy_to_mat('189','193','/data/URY/segy/2013_pgs_uruguay_processing/BG_Total_Merge_Working/PreMatch_Acors2.segy');

%set some constants
start_samp = 1;
end_samp = 251;
max_lag = 125;
weiner = 1;
showqc = 0;

if weiner == 1
    %==========================================================================
    % test of a weiner filter applied to a zero phase wavelet made from the
    % freq spectra of the auto-correlations
    % take auto correlations and make wavelets
    % constants ===============================================================
    NF = 91; %matching filter length
    percenttaper = 20;  % taper to apply to the start and end of wavelet
    % e.g 5 is 5% of the length of the wavelet at the start and at the end
    percenttaperau = 20;
    %==========================================================================
    %tapered_acors = bsxfun(@times,input_acors,taper);
    figure(20); plot(input_acors);
    tapered_acors = taperwavelet(input_acors,percenttaperau);
    figure(21); plot(tapered_acors);
    
    nfft=2^nextpow2(size(tapered_acors,1)); % find next power of 2 for the fft
    
    input_ffts = fft(tapered_acors,nfft); % fft the autocorrelations
    
    freq_samps=nfft/2; % fft outputs a "mirrored" transform so we only need to keep half of it
    
    input_specs = sqrt(abs(input_ffts(1:freq_samps,:))); % amp spec is complex magnitude of transform
    
    % average the amplitude spectra
    %input_specs_avg = mean(input_specs(:,1:2:end),2);
    %input_specs_avg = [input_specs_avg mean(input_specs(:,2:2:end),2)];
    
    %figure(22); plot(input_specs_avg);
    %filt_specs = [input_specs_avg;flipud(input_specs_avg)]; % flip and duplicate spectrum for ifft
    filt_specs = [input_specs;flipud(input_specs)]; % flip and duplicate spectrum for ifft
    
    zerowavelets = ifft(filt_specs,'symmetric'); % transform back to time domain
    zerowavelets=[zerowavelets(freq_samps+1:freq_samps*2,:);zerowavelets(1:freq_samps,:)]; % fix the time zero of the filter
    
    % now take the wavelets and make weiner filter to match them
    wref = zerowavelets(:,1:2:end);
    wout = zerowavelets(:,2:2:end);
    
    wref = taperwavelet(wref,percenttaper);
    wout = taperwavelet(wout,percenttaper);
    
    NF = (length(wref)*2)-1;
    
    filts = zeros(NF,size(wref,2));
    output = zeros(size(wref,1),size(wref,2));;
    
    for fl = 1:size(wref,2)
        %   [a,b] = calculate_filter(wref(:,fl),wout(:,fl),NF,1);
        [filts(:,fl),output(:,fl)] = calculate_filter(wref(:,fl),wout(:,fl),NF,1);
    end
    
    %avfftfilt = amp_spec(filter,4);
    %avfftwref = amp_spec(wref,4);
    %avfftoutput = amp_spec(output,4);
    
    
    if showqc == 1
        figure(2)
        axis tight
        subplot(3,2,1); plot(wref)
        title('input wavelet')
        subplot(3,2,2); plot(wout)
        title('desired wavelet')
        subplot(3,2,3); plot(filter)
        title('weiner filter')
        subplot(3,2,4); %plot(wref,'-b');
        title('input wavelet in blue, desired in green, matched in red')
        plot(wout,'-g');
        hold all;  plot(output,'--r');
        hold off
        subplot(3,2,5); amp_spec(filter,4);
        xlabel('Frequency (Hz)')
        ylabel('Power dB')
        title('Weiner filter amp spectrum')
        subplot(3,2,6); amp_spec(wout,4);
        hold all;
        amp_spec(output,4);
        hold off
        xlabel('Frequency (Hz)')
        ylabel('Power dB')
        title('Amplitude spectrum, desired blue, output green')
    end
else
    %==========================================================================
    % this part is using spectral matching
    taper = 0.5*(1+cos(-pi/max_lag.*[-max_lag:max_lag]))'; % calculate cosine taper function
    
    tapered_acors = bsxfun(@times,input_acors,taper);
    
    figure(23); plot(tapered_acors);
    nfft=2^nextpow2(size(tapered_acors,1)); % find next power of 2 for the fft
    freq_srate=(500/(input_str.s_rate/1000))/(nfft/2); % find the frequency sample interval from Nyquist and number of samples
    
    input_ffts = fft(tapered_acors,nfft); % fft the autocorrelations
    
    freq_samps=nfft/2; % fft outputs a "mirrored" transform so we only need to keep half of it
    
    input_specs = sqrt(abs(input_ffts(1:freq_samps,:))); % amp spec is complex magnitude of transform
    
    
    
    
    input_specs_pw = 0.01*mean(input_specs); % add prewhitening to stabilise filter at frequencies with v low amplitude
    
    survey1_specs = bsxfun(@plus,input_specs(:,1:2:9),input_specs_pw(1:2:9));
    survey2_specs = bsxfun(@plus,input_specs(:,2:2:10),input_specs_pw(1:2:9));
    
    filt_specs = survey2_specs./survey1_specs; % divide desired spectrum by input to get filter spectrum
    
    filt_specs = [filt_specs;flipud(filt_specs)]; % flip and duplicate spectrum for ifft
    
    filts = ifft(filt_specs,'symmetric'); % transform back to time domain
    filts=[filts(freq_samps+1:freq_samps*2,:);filts(1:freq_samps,:)]; % fix the time zero of the filter
    
    figure(26); plot(filt_specs);
    figure(27); plot(filts);
    % set up output
end

output_vol{1,1} = 'Meta data for output files';
output_vol{2,2} = filts; % store output traces in inline/xline order
output_vol{1,2}{1,1} = ilxl(1:2:9,:); % store inline/xline values
output_vol{1,2}{2,1} = uint32(zeros(5)); % what is this for? offsets?
output_vol{1,1} = strcat('amp matching filters ',date);
output_vol{1,3} = 'is_gather'; % 1 is yes, 0 is no
output_vol{2,3}=0;
output_vol{2,1} = 'amp_match_filters_wf_cj';

output_dir = '/data/URY/segy/2013_pgs_uruguay_processing/BG_Total_Merge_Working/';
% check to see if the directory exists
if exist(output_dir,'dir') == 0
    mkdir(output_dir);
end



node_segy_write(output_vol,'1',input_str.s_rate/1000,output_dir); % write the output

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

