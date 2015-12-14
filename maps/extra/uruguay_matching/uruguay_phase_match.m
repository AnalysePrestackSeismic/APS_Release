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
% 2. Read in average crosscorrelations
% 3. Taper the average crosscorrelation with a cosine taper
% 4. Force amp spec to 1 by making complex magnitude equal to 1
% 7. Inverse fft to get filter 

[input_str ilxl input_xcors] = segy_to_mat('189','193','/data/URY/segy/2013_pgs_uruguay_processing/BG_Total_Merge_Working/PreMatch_Xcors2.segy');


start_samp = 1;
end_samp = 101;
max_lag = 50;

taper = 0.5*(1+cos(-pi/max_lag.*[-max_lag:max_lag]))'; % calculate cosine taper function

tapered_xcors = bsxfun(@times,input_xcors,taper);



nfft=2^nextpow2(size(tapered_xcors,1)); % find next power of 2 for the fft
freq_srate=(500/(input_str.s_rate/1000))/(nfft/2); % find the frequency sample interval from Nyquist and number of samples

input_ffts = fft(tapered_xcors,nfft); % fft the autocorrelations

freq_samps=nfft/2; % fft outputs a "mirrored" transform so we only need to keep half of it

amp_specs = ones(size(input_ffts)); % set amp spec of filters to 1
phase_specs = angle(input_ffts);
phase_specs = -phase_specs; % in this case, reverse sign of phase spec because want to match to second trace in xcor

filt_ffts = amp_specs.*exp(1i*phase_specs); % convert amp and phase spec back to complex series


filts = ifft(filt_ffts,'symmetric'); % transform back to time domain
filts=[filts(freq_samps+1:freq_samps*2,:);filts(1:freq_samps,:)]; % fix the time zero of the filter 


% set up output

output_vol{1,1} = 'Meta data for output files';
output_vol{2,2} = filts; % store output traces in inline/xline order
output_vol{1,2}{1,1} = ilxl(1:2:9,:); % store inline/xline values
output_vol{1,2}{2,1} = uint32(zeros(5)); % what is this for? offsets?
output_vol{1,1} = strcat('phase matching filters ',date);
output_vol{1,3} = 'is_gather'; % 1 is yes, 0 is no
output_vol{2,3}=0;
output_vol{2,1} = 'phase_match_filters';

output_dir = '/data/URY/segy/2013_pgs_uruguay_processing/BG_Total_Merge_Working/';
% check to see if the directory exists
if exist(output_dir,'dir') == 0
    mkdir(output_dir);
end

figure; plot(filt_specs);
figure; plot(filts);

node_segy_write(output_vol,'1',input_str.s_rate/1000,output_dir); % write the output



