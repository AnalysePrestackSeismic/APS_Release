function [ peak_freqs ] = peak_freq_calc( job_meta,traces )
% Pick dominant frequency
%  
% Input traces would normally a stack volume
%
% Take FFTs of input traces over small windows
% Smooth the amplitude spectra
% Select frequencies where amplitude > 0.75 * maximum
% 
% peak_freqs has 4 elements:
% - time of window centre
% - lowest frequency with amplitude over 75% of max
% - highest frequency with amplitude over 75% of max
% - weighted avg of lowest and highest (weighted to highest)


% Hardwired parameters:

ns_win = 128;
fft_samps = ns_win/2;
ns_overlap = 96;
taperlen = 16;
filt_smof =  ones(1,9)/9;
freq_thresh = 0.65;

% calculate some more parameters:

[n_samples,~] = size(traces);

start_index = 1:ns_win-ns_overlap-1:n_samples-ns_win;
end_index = start_index+ns_win-1;
window_centre = start_index+fft_samps;

freq_axis = (1e6/job_meta.s_rate)/2*linspace(0,1,ns_win/2);

num_windows = length(start_index);
peak_freqs = zeros(num_windows,4);

taperst = (sin(linspace((-pi/2),(pi/2),taperlen)')+1)/2;
taperend = 1 - taperst;
taperapply = [taperst;ones((ns_win-(taperlen*2)),1);taperend];

% calculate peak frequencies 

for ii = 1:num_windows
    fft_out = zeros(2,2);
    fft_out = abs(fft(bsxfun(@times,traces(start_index(ii):end_index(ii),:),taperapply)));
    avgfreq = sum(fft_out,2,'double');
    avgfreqsmo = conv(avgfreq(1:fft_samps),filt_smof,'same');
    high_amp_mask = (avgfreqsmo > max(avgfreqsmo)*0.75);
    [~,rhs] = max(high_amp_mask(end:-1:1));
    [~,lowf_idx] = max(high_amp_mask(1:1:end));
    highf_idx = (fft_samps - rhs) + 1;
   
    peak_freqs(ii,1:3) = [window_centre(ii) freq_axis(lowf_idx) freq_axis(highf_idx)];
    peak_freqs(ii,4) = (1-freq_thresh)*peak_freqs(ii,2) + freq_thresh*peak_freqs(ii,3);
    
    end

end

