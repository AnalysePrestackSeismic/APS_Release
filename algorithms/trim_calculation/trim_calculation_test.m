function [] = trim_calculation_test(job_meta_path,i_block,startvol,endvol,shift,zsmooth,itm,otm)
%function [] = trim_calculation_test(job_meta_path,i_block,startvol,endvol,tottracerun,maxzout,shift,zsmooth)
%
% this works out and applies trim statics and keeps a sum of the statics
% for each z sample
% ============================================================================================================
% sort out parameters
startvol = str2double(startvol);
endvol = str2double(endvol);
otm = str2double(otm);
itm = str2double(itm);
shift = str2double(shift);
zsmooth = str2double(zsmooth);
zsmooth2 = shift*2;
if zsmooth2 < 8
    zsmooth2 = 8;
end    
shiftinc = 0.5;
% ============================================================================================================
% load the job meta file containing more paramters
job_meta = load(job_meta_path);
sample_rate = (job_meta.s_rate/1000);
% add the history of jobs run and this one to the curent ebcdic
ebdichdr = ['trim statics '];
if isfield(job_meta,'comm_history')
    ebdichdr2 = job_meta.comm_history;
    tmpebc = ebdichdr2{size(ebdichdr2,1),2};
else
    ebdichdr2{1,2} = '';
    tmpebc = '';
end

for ebcii = (size(ebdichdr2,1)-1):-1:1
    tmpebcc = regexp(ebdichdr2{ebcii,2},'/','split');
    tmpebc = [tmpebc tmpebcc{1}  tmpebcc{end}];
end
tmpebc = sprintf('%-3200.3200s',tmpebc);
clear tmpebcc ebdichdr2;
%
% ============================================================================================================
% Read the data for this block

% [~, traces, ilxl_read, offset_read] = node_segy_read(job_meta_path,'1',i_block);
% offset = unique(offset_read);
% traces = reshape(traces,size(traces,1),length(offset),[]);
% 
% [n_samples,n_traces_gather,n_traces] = size(traces);

[~, results_out{3,2}, results_out{1,2}{1,1}, offset_read] = node_segy_read(job_meta_path,'1',i_block);
offset = unique(offset_read);
results_out{3,2} = reshape(results_out{3,2},size(results_out{3,2},1),length(offset),[]);

% drop data as required
% order is samples, angles, skeys
% results_out{3,2} = results_out{3,2}(1:maxzout,:,:);
% job_meta.n_samples{1} = maxzout;

%results_out{1,2}{1,1}

[n_samples,n_traces_gather,n_traces] = size(results_out{3,2});

% apply any itm and otm to the gathers
if itm > 0
    results_out{3,2}(:,1:itm,:) = 0;
end
if otm > 0
    results_out{3,2}(:,otm:end,:) = 0;
end

% ============================================================================================================
% Make results meta information

% for the pre-stack dataset
results_out{1,1} = 'Meta data for output files';
%results_out{resultno,2}{1,1} = ilxl_read;
results_out{1,2}{2,1} = offset_read';
ebcstrtowrite = sprintf('%-3200.3200s',[results_out{1,1} '  ' ebdichdr '  ' tmpebc]);
results_out{1,1} = ebcstrtowrite;

% for the stack output
results_out2{1,1} = 'Meta data for output files';
results_out2{1,2}{1,1} = results_out{1,2}{1,1}(1:n_traces_gather:end,:);
results_out2{1,2}{2,1} = uint32(zeros(n_traces,1));
%ebcstrtowrite = sprintf('%-3200.3200s',[results_out{resultno,1} '  ' ebdichdr '  ' tmpebc]);
results_out2{1,1} = ebcstrtowrite;

output_dir = [job_meta.output_dir,'trimout/'];
% prestack output
filename = 'gaths';
results_out{2,1} = strcat(filename,'_trim_shifts_',num2str(shift),'-',num2str(zsmooth));
results_out{3,1} = strcat(filename,'_trim_data_',num2str(shift),'-',num2str(zsmooth));

results_out{2,2} = zeros(n_samples,n_traces_gather,n_traces,'int16');
%results_out{3,2} = zeros(n_samples,n_traces_gather,n_traces,'single');

% stack result
filename2 = 'stack';
results_out2{2,1} = strcat(filename2,'_trim_sum_',num2str(shift),'-',num2str(zsmooth));
results_out2{2,2} = zeros(n_samples,n_traces,'single');
results_out2{3,1} = strcat(filename2,'_pretrim_',num2str(shift),'-',num2str(zsmooth));
results_out2{3,2} = zeros(n_samples,n_traces,'single');
results_out2{4,1} = strcat(filename2,'_posttrim_',num2str(shift),'-',num2str(zsmooth));
results_out2{4,2} = zeros(n_samples,n_traces,'single');


% ============================================================================================================
% loop round the traces 3d array one gather at a time and do calculations
f_max = (1/sample_rate)*1000;
time = repmat((0:sample_rate:(n_samples-1)*sample_rate)',1,n_traces_gather);
phase_shifts = bsxfun(@times,(-shift:shiftinc:shift),(1/1000).*2.*pi.*repmat((0:f_max/(n_samples-1):f_max)',1,1+2*(shift/shiftinc)));
totnumshifts = length(-shift:shiftinc:shift);
S = (1/zsmooth)*spdiags(repmat([(1:1:zsmooth),(zsmooth-1:-1:1)],n_samples,1),[(-zsmooth+1:1:0),(1:1:zsmooth-1)],n_samples,n_samples);
S2 = (1/zsmooth2)*spdiags(repmat([(1:1:zsmooth2),(zsmooth2-1:-1:1)],n_samples,1),[(-zsmooth2+1:1:0),(1:1:zsmooth2-1)],n_samples,n_samples);

%filttraces = [1 2 3 2 1]/9;
filt_smo =  ones(1,5)/5;

for kk = 1:n_traces;

    % initialise zeros in output
    trim_shift = zeros(n_samples,n_traces_gather,'single');
    trim_shiftb = zeros(n_samples,n_traces_gather,'single');
    
    % get the input gather
    trim_data = results_out{3,2}(:,:,kk);   
    
    % apply an automatic mute to remove low amp noise
    mask = low_amp_mute(trim_data);
    trim_data = trim_data .* mask;
    
    % balence the data to avoid any loud parts dominating any solution
    trim_data_filt = time_balence(trim_data);
    %trim_data_filt = trim_data;  
    
    figure(101); imagesc(trim_data_filt); colormap(gray); caxis([-4000 4000]);
    
    % apply a small smoother to reduce noise
    %trim_data_filt = conv2(1,filt_smo,trim_data_filt,'same');
    trim_data_filt = medfilt3(trim_data_filt,[0 5],0);
    figure(106); imagesc(trim_data_filt); colormap(gray); caxis([-4000 4000]);
    trim_data_filt = trim_data_filt .* mask;
    
    % flip the data to work from outside in
    trim_data_filt = fliplr(trim_data_filt);

    for ii = 2:n_traces_gather
        t1 = double(trim_data_filt(:,ii-1));
        T1 = S*spdiags(t1,0,n_samples,n_samples);
        T1b = S2*spdiags(t1,0,n_samples,n_samples);
        
        t2 = double(trim_data_filt(:,ii));
        t2_shifts = ifft(bsxfun(@times,fft(t2),exp(1i*phase_shifts)),'symmetric');
        
        for count = 1:totnumshifts
            T2 = S*spdiags(t2_shifts(:,count),0,n_samples,n_samples);
            T2b = S2*spdiags(t2_shifts(:,count),0,n_samples,n_samples);
            det_coef(:,count) =  ((T1*t2_shifts(:,count)).*(T2*t1))./((T1*t1).*(T2*t2_shifts(:,count)));
            det_coefb(:,count) =  ((T1b*t2_shifts(:,count)).*(T2b*t1))./((T1b*t1).*(T2b*t2_shifts(:,count)));
        end
        
        %[~,idx] = max(det_coef');
        %trim_shift(:,ii-1) = idx'-shift-1;      
        [~,idx]= max(det_coef,[],2);
        [~,idxb]= max(det_coefb,[],2);
        
        trim_shift(:,ii) = idx-(shift/shiftinc)-1;
        trim_shiftb(:,ii) = idxb-(shift/shiftinc)-1;
        %trim_shift(:,ii-1) = idx-shift-1;
        
    end
    
    maskb = [ones(n_samples,1,'single'),mask(:,1:(end -1))] .* [mask(:,2:end),ones(n_samples,1,'single')];
    trim_shift = [trim_shift(:,1),fliplr(trim_shift(:,2:end))].* maskb;
    trim_shiftb = [trim_shiftb(:,1),fliplr(trim_shiftb(:,2:end))].* maskb;
    %trim_shift = trim_shift .* maskb;
    %trimfold = max(cumsum(maskb,2),[],2);
    %figure;
    figure(100); imagesc(trim_shift); caxis([-5 5]);
    %imagesc(trim_shiftb); 
    
    %trim_sum = sum(abs(trim_shift(:,startvol:endvol)),2) ./ max(cumsum(maskb(:,startvol:endvol),2),[],2);
    
    trim_shift = cumsum(trim_shift,2);
    trim_shift = trim_shift .* mask;
    %trim_shift(data==0) = 0;
    %         trim_sum = sum(trim_shift,2);
    %n = length(trim_shift);
    %trim_sum = sqrt((1/n)*sum((trim_shift.^2),2));
    stackin = sum(trim_data(:,startvol:endvol),2) ./ max(cumsum(maskb(:,startvol:endvol),2),[],2);
    
    for ii = 1:n_traces_gather
        trim_data(:,ii) = interp1(time(:,ii),trim_data(:,ii),time(:,ii)-trim_shift(:,ii),'linear',0);
        trim_shiftb_out(:,ii) = interp1(time(:,ii),trim_shiftb(:,ii),time(:,ii)-trim_shift(:,ii),'linear',0);
    end
    
    figure(102); imagesc(time_balence(trim_data)); colormap(gray); caxis([-4000 4000]);   
    figure(103); imagesc(trim_data); colormap(gray);
    
    % normalise the sum of the shifts for the fold in case we have
    % difference in fold
    trim_sum = sum(abs(trim_shiftb_out(:,startvol:endvol)),2) ./ max(cumsum(maskb(:,startvol:endvol),2),[],2);
    %trim_sum = sum(abs(trim_shiftb_out(:,startvol:endvol)),2);
    figure(109); imagesc(trim_shiftb_out); caxis([-5 5]);
    
    stackout = sum(trim_data(:,startvol:endvol),2) ./ max(cumsum(maskb(:,startvol:endvol),2),[],2);
    
    results_out{2,2}(:,:,kk) = int16(trim_shift);
    results_out{3,2}(:,:,kk) = trim_data;
    results_out2{2,2}(:,kk) = trim_sum;
    results_out2{3,2}(:,kk) = stackin;
    results_out2{4,2}(:,kk) = stackout;
    
end

% need to reshape the 2 3d matricies into 2d ones as the segy write just wants samples * total traces

results_out{2,2} = reshape(results_out{2,2},n_samples,[]);
results_out{3,2} = reshape(results_out{3,2},n_samples,[]);

% check to see if the directory exists
if exist(output_dir,'dir') == 0
    mkdir(output_dir);    
end

i_block = str2double(i_block);
% write the pre-stack dataset
node_segy_write(results_out,i_block, sample_rate, output_dir)
% write the stack
node_segy_write(results_out2,i_block, sample_rate, output_dir)

end