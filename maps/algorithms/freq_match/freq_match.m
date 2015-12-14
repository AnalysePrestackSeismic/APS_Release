function [  ] = freq_match( ref_stack_meta_path, i_block, input_stack_meta_path, output_path, start_win, end_win, max_lag )
%Design and apply frequency-matching filter
%   
%Measure freqency content of a reference stack and input stack and design a
%zero-phase filter to match them
%
% ref_stack_meta_path:      meta data for reference stack
% i_block:                  block number to analyse (0 = analyse all blocks)
% input_stack_meta_path:    meta data for input stack
% output path:              path for output segy dataset
% start_win:                start window
% end_win:                  end window
% max_lag:                      max autocorrelation lead/lag

warning off all;

ref_meta = load(ref_stack_meta_path);
input_meta = load(input_stack_meta_path);

ref_acor_sum = zeros(max_lag*2+1,1);

start_samp = round(start_win/(ref_meta.s_rate/1000));
end_samp = round(end_win/(ref_meta.s_rate/1000));
win_len = end_samp - start_samp + 1;
    
for ii=1:size(ref_meta.liveblocks,1)
    
    i_block=int2str(ref_meta.liveblocks(ii));
    
    [ref_seismic, ref_traces, ~, ~] = node_segy_read(ref_stack_meta_path,'1',i_block);
    ref_ntraces=size(ref_traces,2);
    ref_acors = zeros(2*max_lag+1,ref_ntraces);
    
    for trace=1:ref_ntraces
        
        ref_acors(:,trace)=xcorr(ref_traces(start_samp:end_samp,trace),max_lag);
        
    end
    
    ref_acor_sum = ref_acor_sum+sum(ref_acors,2);
    
end

input_acor_sum = zeros(max_lag*2+1,1);

for ii=1:size(input_meta.liveblocks,1)
    
    i_block=int2str(input_meta.liveblocks(ii));
    
    
    [input_seismic, input_traces, ~, ~] = node_segy_read(input_stack_meta_path,'1',i_block);
    input_ntraces=size(input_traces,2);
    input_acors = zeros(2*max_lag+1,input_ntraces);
    
    
    for trace=1:input_ntraces
        
        input_acors(:,trace)=xcorr(input_traces(start_samp:end_samp,trace),max_lag);
        
    end
    
    
    input_acor_sum = input_acor_sum + sum(input_acors,2);
    
end


taper = 0.5*(1+cos(-pi/max_lag.*[-max_lag:max_lag]))';
ref_acor_sum = ref_acor_sum.*taper;
input_acor_sum = input_acor_sum.*taper;

nfft=2^nextpow2(size(ref_acor_sum,1));
freq_srate=(500/(ref_seismic.s_rate/1000))/(nfft/2);

ref_fft = fft(ref_acor_sum,nfft);
input_fft = fft(input_acor_sum,nfft);

freq_samps=nfft/2;

ref_spec = sqrt(abs(ref_fft(1:freq_samps)));
input_spec = sqrt(abs(input_fft(1:freq_samps)));

ref_pw = 0.01*mean(ref_spec);
input_pw = 0.01*mean(input_spec);

filt_spec = (ref_spec+ref_pw)./(input_spec+input_pw);
filt_spec = filt_spec./mean(filt_spec);

filt_spec = [filt_spec;flipud(filt_spec)];

filt = ifft(filt_spec,'symmetric');
filt=[filt(ceil(freq_samps/2):freq_samps);filt(1:floor(freq_samps/2))];


output_traces=zeros(input_seismic.n_samples,input_seismic.n_traces);
output_ilxl=zeros(input_seismic.n_traces,3);
count=1;

for ii=1:size(input_meta.liveblocks,1)
    
    i_block=int2str(input_meta.liveblocks(ii));
    
    
    [input_seismic, input_traces, input_ilxl, ~] = node_segy_read(input_stack_meta_path,'1',i_block);
    input_ntraces=size(input_traces,2);
    output_traces(:,count:count+input_ntraces-1) = convn(input_traces,filt,'same');
    output_ilxl(count:count+input_ntraces-1,1) = count:count+input_ntraces-1;
    output_ilxl(count:count+input_ntraces-1,2:3) = input_ilxl;
    count=count+input_ntraces;
    
end

sort_order = sortrows(output_ilxl,[2 3]);
sort_order = sort_order(:,1);
sort_order = sort_order(sort_order>0)';

output_vol{1,1} = 'Meta data for output files';
output_vol{2,2} = output_traces(:,sort_order);
output_vol{1,2}{1,1} = output_ilxl(sort_order',2:3);
output_vol{1,2}{2,1} = uint32(zeros(size(sort_order,2),1)); % what is this for? offsets?
output_vol{1,1} = strcat('Output from freq matching ',date);
output_vol{1,3} = 'is_gather'; % 1 is yes, 0 is no
output_vol{2,3}=0;
output_vol{2,1} = 'freq_match_stack';

output_dir = [input_meta.output_dir,'freq_match/'];
% check to see if the directory exists
if exist(output_dir,'dir') == 0
    mkdir(output_dir);
end


% save('/segy/BOL/LVT_Los_Suris_2014/matlab/freq_match_test.mat');



node_segy_write(output_vol,'1',input_seismic.s_rate/1000,output_path);


end

