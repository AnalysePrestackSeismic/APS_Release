function [  ] = freq_scale( stack_meta_path, i_block, output_path, start_win, end_win, freq_pairs, max_lag )
%Design and apply frequency-matching filter
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
% 1. Autocorrelate every trace in input using specified window and max lag
% 2. Average the autocorrelations
% 3. Taper the average autocorrelation with a cosine taper
% 4. Take square root of amplitude spectrum
% 5. Interpolate supplied freq,amp pairs onto same freq sampling as above
% 6. Divide output spectrum by input spectrum to get filter spectrum
% 7. Inverse fft to get filter response
% 8. Convolve with input traces and write output to segy

warning off all;

stack_meta = load(stack_meta_path); % load meta data for input stack

start_samp = round(start_win/(stack_meta.s_rate/1000)); % convert window start to samples
end_samp = round(end_win/(stack_meta.s_rate/1000)); % convert window end to samples
    
stack_acor_sum = zeros(max_lag*2+1,1); % initialise matrix for average autocorrelation
total_traces = 0; % counter for total number of traces

for ii=1:size(stack_meta.liveblocks,1) % loop round the live blocks in the input volume
    
    i_block=int2str(stack_meta.liveblocks(ii)); % convert i_block to string for node_segy_read
        
    [stack_seismic, stack_traces, ~, ~] = node_segy_read(stack_meta_path,'1',i_block); % read seismic and meta data
    stack_ntraces=size(stack_traces,2); % find number of traces in the block
    total_traces = total_traces + stack_ntraces; % update total number of traces
            
    for trace=1:stack_ntraces % loop round the traces in the block
        % accumulate each autocorrelation into the sum
        stack_acor_sum = stack_acor_sum + xcorr(stack_traces(start_samp:end_samp,trace),max_lag); 
        
    end
           
end

taper = 0.5*(1+cos(-pi/max_lag.*[-max_lag:max_lag]))'; % calculate cosine taper function
stack_acor_sum = (stack_acor_sum.*taper)./total_traces; % apply taper and divide by number of traces

nfft=2^nextpow2(size(stack_acor_sum,1)); % find next power of 2 for the fft
freq_srate=(500/(stack_seismic.s_rate/1000))/(nfft/2); % find the frequency sample interval from Nyquist and number of samples

stack_fft = fft(stack_acor_sum,nfft); % fft the autocorrelation

freq_samps=nfft/2; % fft outputs a "mirrored" transform so we only need to keep half of it

stack_spec = sqrt(abs(stack_fft(1:freq_samps))); % amp spec is complex magnitude of transform

% map supplied freq values onto freq samples

freq_samp_pairs(:,1) = freq_pairs(:,1)./freq_srate;
freq_samp_pairs(:,2) = freq_pairs(:,2);

% interpolate spectrum for all frequencies

ref_spec = interp1([1;freq_samp_pairs(:,1);freq_samps],[0;freq_samp_pairs(:,2);0],[1:freq_samps]');

ref_pw = 0.01*mean(ref_spec); % add prewhitening to stabilise filter at frequencies with v low amplitude
stack_pw = 0.01*mean(stack_spec);

filt_spec = (ref_spec+ref_pw)./(stack_spec+stack_pw); % divide desired spectrum by input to get filter spectrum
filt_spec = filt_spec./mean(filt_spec); % normalise the filter

filt_spec = [filt_spec;flipud(filt_spec)]; % flip and duplicate spectrum for ifft

filt = ifft(filt_spec,'symmetric'); % transform back to time domain
filt=[filt(freq_samps+1:freq_samps*2);filt(1:freq_samps)]; % fix the time zero of the filter 

% now apply the filter to the input traces
% also need to sort the traces into inline/xline order
% this is done by making an array of the input inline/xline values with an
% extra column for the index for each trace
% then use sortrows to sort into inline/xline and the other column then
% contains a vector you can use to index into the trace array

output_traces=zeros(stack_seismic.n_samples,stack_seismic.n_traces); % initialise matrix for output
output_ilxl=zeros(stack_seismic.n_traces,3); % initialise counter/inline/xline array for output
count=1; % counter to keep track of where we are in the output volume

for ii=1:size(stack_meta.liveblocks,1) % loop round the live blocks
    
    i_block=int2str(stack_meta.liveblocks(ii)); % convert i_block to string for node_segy_read
        
    [stack_seismic, stack_traces, stack_ilxl, ~] = node_segy_read(stack_meta_path,'1',i_block); % read the data
    stack_ntraces=size(stack_traces,2); % calculate number of traces in the block
    output_traces(:,count:count+stack_ntraces-1) = convn(stack_traces,filt,'same'); % apply the filter
    output_ilxl(count:count+stack_ntraces-1,1) = count:count+stack_ntraces-1; % store the counter for sorting later 
    output_ilxl(count:count+stack_ntraces-1,2:3) = stack_ilxl; % store inline/xline values
    count=count+stack_ntraces; % update the counter
    
end

% output_acor_sum = zeros(size(stack_acor_sum));
% 
% for trace=1:(count-1)
%     output_acor_sum=output_acor_sum+xcorr(output_traces(start_samp:end_samp,trace),max_lag);
% end

sort_order = sortrows(output_ilxl,[2 3]); % sort the inline/xline array
sort_order = sort_order(:,1); % make an array of just the counter to use as an index
sort_order = sort_order(sort_order>0)'; % discard any zero values

% set up output

output_vol{1,1} = 'Meta data for output files';
output_vol{2,2} = output_traces(:,sort_order); % store output traces in inline/xline order
output_vol{1,2}{1,1} = output_ilxl(sort_order',2:3); % store inline/xline values
output_vol{1,2}{2,1} = uint32(zeros(size(sort_order,2),1)); % what is this for? offsets?
output_vol{1,1} = strcat('Output from freq matching ',date);
output_vol{1,3} = 'is_gather'; % 1 is yes, 0 is no
output_vol{2,3}=0;
output_vol{2,1} = 'freq_scale_stack';

output_dir = [stack_meta.output_dir,'freq_scale/'];
% check to see if the directory exists
if exist(output_dir,'dir') == 0
    mkdir(output_dir);
end



node_segy_write(output_vol,'1',stack_seismic.s_rate/1000,output_path); % write the output


end

