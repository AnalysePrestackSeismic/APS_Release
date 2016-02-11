function [out_traces stack_shifts stack_traces coeff_sum ] = trim_shifts_calculate(mode,srate,traces,timbal,shift_param)
% this works out and applies trim statics and keeps a sum of the statics
% for each z sample
% 
% mode: 'calc' or 'apply'
% job_meta: meta data for the input gathers
% traces: input traces - should be cdp-sorted nmo-corrected muted gathers
% timbal: flag for time-balancing
% zsmooth: time/depth smoothing
% 
% out_traces: traces with shifts applied or shifts
% out_shifts: traces with samples replaced by shift values
% xcor_coeffs: one trace per input gather with average xcor coefficient
% 
% ============================================================================================================
% sort out parameters
timbal = str2double(timbal);
shift_param = str2double(shift_param);
% ============================================================================================================
% zsmooth parameter was not used

% useselectemode = 0;
% if isempty(regexp(zsmooth,'il','once')) == 0
%     useselectemode = 1;
%     requiredinline =  str2double(regexprep(zsmooth,'il',''));
%     %requiredinline =  str2double(strrep(tottracerun,'il',''));
%     zsmooth = 20;
% else
%     zsmooth = str2double(zsmooth);
% end

% zsmooth = str2double(zsmooth);
%zsmooth2 = shift*2;
% if zsmooth2 < 7
%     zsmooth2 = 7;
% end
% ============================================================================================================


shiftinc = 0.2;
plots = 0;
padding = 50;
% Wavelet estimation parameters
ns_win = 128;
hns_win = ns_win/2;
ns_overlap = 96;
totmaxshift = 5;

% ============================================================================================================
% % load the job meta file containing more parameters
% job_meta = load(job_meta_path);
%sample_rate = (job_meta.s_rate/1000);
sample_rate = srate/1000;

% add the history of jobs run and this one to the curent ebcdic
% ebdichdr = ['trim statics '];
% if isfield(job_meta,'comm_history')
%     ebdichdr2 = job_meta.comm_history;
%     tmpebc = ebdichdr2{size(ebdichdr2,1),2};
% else
%     ebdichdr2{1,2} = '';
%     tmpebc = '';
% end
% 
% for ebcii = (size(ebdichdr2,1)-1):-1:1
%     tmpebcc = regexp(ebdichdr2{ebcii,2},'/','split');
%     tmpebc = [tmpebc tmpebcc{1}  tmpebcc{end}];
% end
% tmpebc = sprintf('%-3200.3200s',tmpebc);
% clear tmpebcc ebdichdr2;
% %
% ============================================================================================================
% Read the data for this block

% [~, traces, ilxl_read, offset_read] = node_segy_read(job_meta_path,'1',i_block);
% offset = unique(offset_read);
% traces = reshape(traces,size(traces,1),length(offset),[]);
%
% [n_samples,n_traces_gather,n_traces] = size(traces);

% [seismic, results_out{3,2}, results_out{1,2}{1,1}, offset_read] = node_segy_read(job_meta_path,'1',i_block);
% offset = unique(offset_read);
% results_out{3,2} = reshape(results_out{3,2},size(results_out{3,2},1),length(offset),[]);
% 
% % check to make sure it read something if not exit
% if isempty(results_out{3,2}) == 1 &&  isempty(results_out{1,2}{1,1}) == 1 && isempty(offset_read) == 1
%     return
% end

% drop data as required
% order is samples, angles, skeys
% results_out{3,2} = results_out{3,2}(1:maxzout,:,:);
% job_meta.n_samples{1} = maxzout;

%results_out{1,2}{1,1}

[n_samples,n_traces_gather,n_traces] = size(traces);
in_n_samples = n_samples;

% apply any itm and otm to the gathers
% sttrc = 2;
% edtrc = n_traces_gather;
% if itm > 0
%     results_out{3,2}(:,1:itm,:) = 0;
%     edtrc = (n_traces_gather - itm);
% end
% if otm > 0
%     results_out{3,2}(:,otm:end,:) = 0;
%     sttrc = (n_traces_gather - otm) + 2;
% end
% midtrc = startvol + floor((endvol - startvol)/2);
% ============================================================================================================
% % Make results meta information
% 
% % for the pre-stack dataset
% results_out{1,1} = 'Meta data for output files';
% %results_out{resultno,2}{1,1} = ilxl_read;
% results_out{1,2}{2,1} = offset_read';
% ebcstrtowrite = sprintf('%-3200.3200s',[results_out{1,1} '  ' ebdichdr '  ' tmpebc]);
% results_out{1,1} = ebcstrtowrite;
% results_out{1,3} = 'is_gather'; % 1 is yes, 0 is no
% 
% 
% % for the stack output
% results_out2{1,1} = 'Meta data for output files';
% results_out2{1,2}{1,1} = results_out{1,2}{1,1}(1:n_traces_gather:end,:);
% results_out2{1,2}{2,1} = uint32(zeros(n_traces,1));
% %ebcstrtowrite = sprintf('%-3200.3200s',[results_out{resultno,1} '  ' ebdichdr '  ' tmpebc]);
% results_out2{1,1} = ebcstrtowrite;
% results_out2{1,3} = 'is_gather'; % 1 is yes, 0 is no
% 
% output_dir = [job_meta.output_dir,'trimout/'];
% % check to see if the directory exists
% if exist(output_dir,'dir') == 0
%     mkdir(output_dir);
% end
% 
% 
% % prestack output
% filename = 'gaths';
% results_out{2,1} = strcat(filename,'_trim_shifts_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
% results_out{3,1} = strcat(filename,'_trim_data_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
% results_out{2,3} = 1;
% results_out{3,3} = 1;
% results_out{2,2} = zeros(in_n_samples,n_traces_gather,n_traces,'single');
% %results_out{3,2} = zeros(n_samples,n_traces_gather,n_traces,'single');
% 
% 
% 
% % stack result
% filename2 = 'stack';
% results_out2{2,1} = strcat(filename2,'_trim_sum_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
% results_out2{2,2} = zeros(in_n_samples,n_traces,'single');
% results_out2{3,1} = strcat(filename2,'_pretrim_guide_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
% results_out2{3,2} = zeros(in_n_samples,n_traces,'single');
% %results_out2{4,1} = strcat(filename2,'_posttrim_no_resid_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
% results_out2{4,1} = strcat(filename2,'_posttrim_guide_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
% results_out2{4,2} = zeros(in_n_samples,n_traces,'single');
% results_out2{5,1} = strcat(filename2,'_posttrim_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
% results_out2{5,2} = zeros(in_n_samples,n_traces,'single');
% results_out2{6,1} = strcat(filename2,'_pretrim_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
% results_out2{6,2} = zeros(in_n_samples,n_traces,'single');
% results_out2{2,3} = 0;
% results_out2{3,3} = 0;
% results_out2{4,3} = 0;
% results_out2{5,3} = 0;
% results_out2{6,3} = 0;
% ========================================================================
% flattern the data to the water bottom and make a stack of the data for
% freq analysis


% % stack to make something to pick wb on
% stackin_freq = sum(results_out{3,2}(:,startvol:midtrc,:),2) ./ ((midtrc - startvol) +1) ;
% % pick wb
% %[wb_idx] = water_bottom_picker(stackin_freq,padding);
% [wb_idx] = water_bottom_picker(squeeze(stackin_freq),padding);
% wb_idx(wb_idx < 0) = 1;
% % write out a wbpick
% %ilxltoprint = results_out{1,2}{1,1}(1:length(offset):end,:);
% % only write out the values that are not 1, so are picked
% %dlmwrite(strcat(output_dir,'wbpick_',i_block,'.xyz'),[ilxltoprint((wb_idx ~= 1)',:),(int32(wb_idx(wb_idx ~= 1)+padding)'.*job_meta.s_rate)/1000],'delimiter', ' ', 'precision', '%-6d','newline', 'unix');
% %
% win_sub = bsxfun(@plus,wb_idx,(0:n_samples-max(wb_idx))');
% win_ind = bsxfun(@plus,win_sub,(0:n_samples:n_samples*(n_traces-1)));
% clear win_sub;
% % make anew stack after flattening, make sure to blank it before
% stackin_freq2 = stackin_freq(win_ind);
% %stackin_freq2 = time_balence_stk(stackin_freq2);
% [n_samples,~] = size(stackin_freq2);
clear win_ind stackin_freq;

wb_idx=ones(n_traces);

stack_fold = squeeze(max(1,sum(abs(traces)>0,2)));

stackin_freq2 = squeeze(sum(traces,2))./(stack_fold/n_traces_gather);

% find the dominant frequency of the data, take a subset of the data and
%
% returns window centre, low, high, weighted avg freq
%
% assign to arrays with same name as in previous code
job_meta.s_rate = srate;

peak_freqs = peak_freq_calc(job_meta,stackin_freq2);

min_freq = peak_freqs(:,[1 2]);
peak_freq = peak_freqs(:,[1 3]);
mid_freq = peak_freqs(:,[1 4]);

peak_freq(:,2) = 90;
mid_freq(:,2) = 80;


% make taper to apply to signal before fft
taperlen = 16;

taperst = (sin(linspace((-pi/2),(pi/2),taperlen)')+1)/2;
taperend = 1 - taperst;
% taperapply = [taperst;ones((ns_win-(taperlen*2)),1);taperend];

% =========================================================================
%see if the code can resample the data down to calculate shifts at
resamp_rt = round(1000 /(max(mid_freq(:,2))*2));
use_samp_drop = 0;
wb_n_samples = n_samples;
orig_sample_rate = sample_rate;
time = repmat((0:sample_rate:(n_samples-1)*sample_rate)',1,n_traces_gather);
if resamp_rt > (sample_rate * 1.8)
    samp_drop = round(resamp_rt/(sample_rate * 1.8) + 0.1);
    sample_rate = sample_rate * samp_drop;
    ufoff = 1000/(sample_rate*2);
    ufon = ufoff*0.8;
    n_samples = length(1:samp_drop:n_samples);
    peak_freq((peak_freq(:,2) > ufoff),2) = ufoff;
    totmaxshift = ceil(sample_rate/2.5); % was 3 before
    use_samp_drop = 1;
end

time_sub = (0:sample_rate:(n_samples-1)*sample_rate);

%
% for jj = 1:gathinc:n_traces;
%     freq_test = results_out{3,2}(:,sttrc:midtrc,jj);
%     %freq_test = time_balence(freq_test);
%
%     for ii = 1:n_win
%         % Estimate wavelets and store meta information
%         % make a tmp array with the fft values in it
%         tmpfft = zeros(2,2);
%         tmpfft = abs(fft(bsxfun(@times,freq_test(start_index(ii):end_index(ii),:),taperapply)));
%         avgfreq = sum(tmpfft,2,'double');
%         avgfreqsmo = conv(avgfreq(1:hns_win),filt_smof,'same');
%         %figure(12); plot(avgfreqsmo);
%         cja = (avgfreqsmo > max(avgfreqsmo)*0.75);
%         [~,rhs] = max(cja(end:-1:1));
%         highf_idx = (hns_win - rhs) + 1;
%         %[~,idxf] = max(avgfreqsmo);
%         %peak_freq(ii,:) = [(start_index(ii)+hns_win) ((peak_freq(ii,2)+freq_axis(highf_idx)))];
%         peak_freq(ii,:) = [(start_index(ii)+hns_win) max(peak_freq(ii,2),freq_axis(highf_idx))];
%     end
% end
%peak_freq(:,2) = peak_freq(:,2)/n_win;
%
% % ========================================================================
% % calculate the shift limits and the smoothing for the trim, currently
% % the wavelength at the peak freq at 2000m/s divided by 30
% % max_shift = [peak_freq(:,1) (2000./(peak_freq(:,2)*42))];

% max_shift = [peak_freq(:,1) (2000./(peak_freq(:,2)*13))];

max_shift = zeros(size(peak_freq));

max_shift(:,1) = peak_freq(:,1);
max_shift(:,2) = shift_param;

% smooth the shifts
%filt_smos =  ones(1,3)/3;
%max_shift = conv(max_shift,filt_smos,'same');
max_shift(:,2) = floor(max_shift(:,2)./shiftinc)*shiftinc;
max_shiftstk = max_shift;
max_shiftstk(:,2) = max_shiftstk(:,2).*4; % was 3 before this allows the residual to do more

% calculate the length of the smoother, currently 10 times the shift value
freq_grid = [max_shift(:,1) floor(max_shift(:,2).*10)];

shift = max(max_shift(:,2));
% find the max shift to calculate
if totmaxshift >= shift
    totmaxshift = shift;
else
    shift = totmaxshift;
end

% clip the max shifts
max_shift((max_shift(:,2) > totmaxshift),2) = totmaxshift;
max_shiftstk((max_shiftstk(:,2) > totmaxshift),2) = totmaxshift;


totnumshifts = length(-shift:shiftinc:shift);

%now make the time varying smoothing as a diagonal on a matrix S
%freq_grid2 = [30 11; 300 21; 900 31];

totalpts = max(freq_grid(:,2));
for mm = 1:size(freq_grid,1)
    freq_interp_grid(mm,:) = [zeros(1,(totalpts-freq_grid(mm,2))) , 0:(totalpts/freq_grid(mm,2)):totalpts ,  (totalpts-(totalpts/freq_grid(mm,2))):-(totalpts/freq_grid(mm,2)):0 , zeros(1,(totalpts-freq_grid(mm,2)))  ];
    shift_interp_grid(mm,:) = [[zeros(1,round(shift/shiftinc-((max_shift(mm,2)/shiftinc)))) , ones(1,round(max_shift(mm,2)/shiftinc))],1,fliplr([zeros(1,round(shift/shiftinc-(max_shift(mm,2)/shiftinc))) , ones(1,round(max_shift(mm,2)/shiftinc))])];
    shift_interp_gridstk(mm,:) = [[zeros(1,round(shift/shiftinc-((max_shiftstk(mm,2)/shiftinc)))) , ones(1,round((max_shiftstk(mm,2))/shiftinc))],1,fliplr([zeros(1,round(shift/shiftinc-((max_shiftstk(mm,2))/shiftinc))) , ones(1,round((max_shiftstk(mm,2))/shiftinc))])];
end
start_interp = 1;
end_interp = n_samples;
freq_zgrid = freq_grid(:,1);

if freq_grid(1,1) > 1
    freq_zgrid = [1; freq_zgrid];
    freq_interp_grid = [ freq_interp_grid(1,:) ; freq_interp_grid];
    shift_interp_grid = [ shift_interp_grid(1,:) ; shift_interp_grid];
    shift_interp_gridstk = [ shift_interp_gridstk(1,:) ; shift_interp_gridstk];
end
if freq_grid(end,1) < n_samples
    freq_zgrid(end+1) = n_samples;
    freq_interp_grid(end+1,:) = freq_interp_grid(end,:);
    shift_interp_grid(end+1,:) = shift_interp_grid(end,:);
    shift_interp_gridstk(end+1,:) = shift_interp_gridstk(end,:);
end

%make the array to go down the diagonal as a smoother
diag_smo_interp = interp1(freq_zgrid,freq_interp_grid,start_interp:1:end_interp,'linear');
S = (1/totalpts)*spdiags(diag_smo_interp,[(-totalpts:1:0),(1:1:totalpts)], n_samples,n_samples);

% now make a smaller smoother
dropf = 3;
dropfsel = [(totalpts + 1) - ((length((totalpts + 1:dropf:(totalpts*2)+1)) -1)*dropf):dropf:totalpts, (totalpts + 1:dropf:(totalpts*2)+1)];
%droprows = [(-round(totalpts/dropf):1:0),(1:1:round(totalpts/dropf))];
droprows = [(-round((length(dropfsel) -1)/2):1:0),(1:1:round((length(dropfsel) -1)/2))];
S2 = (1/totalpts)*spdiags(diag_smo_interp(:,dropfsel),droprows, n_samples,n_samples);

diag_smo_interp3 = interp1(freq_zgrid,freq_interp_grid,start_interp:1:wb_n_samples,'linear');
%S3 = (1/totalpts)*spdiags(diag_smo_interp3(:,dropfsel),droprows, wb_n_samples,wb_n_samples);
S3 = (1/totalpts)*spdiags(diag_smo_interp3,[(-totalpts:1:0),(1:1:totalpts)], wb_n_samples,wb_n_samples);

% make a mask to apply to the shifts before picking the peak
shift_mask = interp1(freq_zgrid,shift_interp_grid,start_interp:1:end_interp,'linear');

shift_interp_gridstk = interp1(2:2:(totnumshifts*2),shift_interp_gridstk',1:((totnumshifts*2)+1),'linear');
shift_interp_gridstk(isnan(shift_interp_gridstk)) = 0;
shift_maskstk = interp1(freq_zgrid,shift_interp_gridstk',start_interp:1:wb_n_samples,'linear');
shift_mask(shift_mask<1) = 0;
shift_maskstk(shift_maskstk<1) = 0;
% ============================================================================================================
% loop round the traces 3d array one gather at a time and do calculations
f_max = (1/sample_rate)*1000;
phase_shifts = bsxfun(@times,(-shift:shiftinc:shift),(1/1000).*2.*pi.*repmat((0:f_max/(n_samples-1):f_max)',1,1+2*(shift/shiftinc)));

f_maxorig = (1/orig_sample_rate)*1000;
phase_shifts_orig = bsxfun(@times,(-shift:shiftinc:shift),(1/1000).*2.*pi.*repmat((0:f_maxorig/(wb_n_samples-1):f_maxorig)',1,1+2*(shift/shiftinc)));

totnumshifts = length(-shift:shiftinc:shift);
%S = (1/zsmooth)*spdiags(repmat([(1:1:zsmooth),(zsmooth-1:-1:1)],n_samples,1),[(-zsmooth+1:1:0),(1:1:zsmooth-1)],n_samples,n_samples);
%S2 = (1/zsmooth2)*spdiags(repmat([(1:1:zsmooth2),(zsmooth2-1:-1:1)],n_samples,1),[(-zsmooth2+1:1:0),(1:1:zsmooth2-1)],n_samples,n_samples);
det_coef = zeros(n_samples,totnumshifts);
det_coefb = zeros(n_samples,totnumshifts);

%filttraces = [1 2 3 3 3 3 3 3 3 2 1]/27;
filttraces = [1 2 3 3 2 1]/12;
filt_smo =  [1 2 1]/4;
filt_smo2 =  ones(1,5)/5;
idxshift = (shift/shiftinc)+1;

% initialise zeros in output
trim_shift = zeros(n_samples,n_traces_gather,'single');
trim_shiftb = zeros(n_samples,n_traces_gather,'single');
trim_shiftb_out = zeros(wb_n_samples,n_traces_gather,'single');
wt_taperapply = [taperst;ones((wb_n_samples-(taperlen*2)),1);taperend];
det_coef3 = zeros(wb_n_samples,length(2:2:totnumshifts) );
norm_fold_b = zeros(wb_n_samples,1,'single');
norm_fold = zeros(wb_n_samples,1,'single'); 
out_traces = zeros(size(traces));
coeff_sum = zeros(wb_n_samples,n_traces);

%
for kk = 1:n_traces;
    %for kk = 1:200;
    %     if useselectemode == 0
    %         requiredinline = results_out{1,2}{1,1}(kk,1);
    %     end
    %     if results_out{1,2}{1,1}(kk,1) == requiredinline;
    
    
    % get the input gather and apply shift to flattern to the water bottom
    trim_data = traces(wb_idx(kk):(wb_idx(kk)+wb_n_samples-1),:,kk);
    
    % apply an automatic mute to remove low amp noise
    fmask = low_amp_mute(trim_data);
    trim_data = trim_data .* fmask;
    
    %drop samples if possible
    if use_samp_drop == 1
        trim_data_filt = bsxfun(@times,trim_data,wt_taperapply);
        trim_data_filt = bandpass_filter(trim_data_filt,(sample_rate/1000),0,0,ufon,ufoff);
        trim_data_filt = trim_data_filt(1:samp_drop:end,:);
        mask = fmask(1:samp_drop:end,:);
    else
        trim_data_filt = trim_data;
        mask = fmask;
    end
    
    % balence the data to avoid any loud parts dominating any solution
    if timbal == 1
        trim_data_filtin = trim_data_filt;
        [trim_data_filt scaleused] = time_balence(trim_data_filt);
    end
    
    %trim_data = trim_data_filt;
    %trim_data_filt = bandpass_filter(d,dt,f1,f2,f3,f4);
    
    if plots == 1; figure(101); imagesc(trim_data_filt); colormap(gray); caxis([-4000 4000]); end;
    
    % apply a small smoother to reduce noise
    %        if itm > 0
    itm=0; otm=n_traces_gather;
    mask(:,1:itm) = 0;
    trim_data_filt(:,1:itm) = repmat(trim_data_filt(:,(itm+1)),1,itm);
    %end
    %if otm > 0
    mask(:,otm:end) = 0;
    trim_data_filt(:,otm:n_traces_gather) = repmat(trim_data_filt(:,(otm-1)),1,(n_traces_gather-otm)+1);
    %end
    
    % apply a small smoother to reduce noise
    %trim_data_filt = medfilt3nt(trim_data_filt,[0 9],0);
%    trim_data_filt = conv2(filt_smo,filttraces,trim_data_filt,'same');
    
    if plots == 1; figure(106); imagesc(trim_data_filt); colormap(gray); caxis([-4000 4000]); end;
    
    %apply mask to the filtered data to remove duplicated data
    trim_data_filt = trim_data_filt .* mask;
    
    % flip the data to work from outside in
    trim_data_filt = fliplr(trim_data_filt);
    
    %for ii = sttrc:edtrc
    for ii = 2:n_traces_gather;
        t1 = double(trim_data_filt(:,ii-1));
        T1 = S*spdiags(t1,0,n_samples,n_samples);
        T1b = S2*spdiags(t1,0,n_samples,n_samples);
        
        %t2 = double(trim_data_filt(:,ii));
        t2f = fft(double(trim_data_filt(:,ii)));
        t2_shifts = ifft(bsxfun(@times,t2f,exp(1i*phase_shifts)),'symmetric');
        
        for count = 1:totnumshifts
            %T2 = S*spdiags(t2_shifts(:,count),0,n_samples,n_samples);
            T2diag = spdiags(t2_shifts(:,count),0,n_samples,n_samples);
            %sparse
            T2 = S*T2diag;
            T2b = S2*spdiags(t2_shifts(:,count),0,n_samples,n_samples);
            det_coef(:,count) =  ((T1*t2_shifts(:,count)).*(T2*t1))./((T1*t1).*(T2*t2_shifts(:,count)));
            det_coefb(:,count) =  ((T1b*t2_shifts(:,count)).*(T2b*t1))./((T1b*t1).*(T2b*t2_shifts(:,count)));
        end
        % add in a taper to zero out large shifts in the shallow before picking
        % the max
        det_coef = det_coef.*shift_mask;
        %det_coefb = det_coefb.*shift_mask;
        % cj edit took out the mask as it was making steps in the
        % results
        
        %[~,idx] = max(det_coef');
        %trim_shift(:,ii-1) = idx'-shift-1;
        [~,idx]= max(det_coef,[],2);
        [coeffb,idxb]= max(det_coefb,[],2);
        
        trim_shift(:,ii) = (idx-idxshift)*shiftinc;
        trim_shiftb(:,ii) = (idxb-idxshift)*shiftinc;
        trim_coeffb(:,ii) = coeffb;
        %trim_shift(:,ii-1) = idx-shift-1;
        
    end
    
    maskb = [ones(n_samples,1,'single'),mask(:,1:(end -1))] .* [mask(:,2:end),ones(n_samples,1,'single')];
    trim_shift = [trim_shift(:,1),fliplr(trim_shift(:,2:end))].* maskb;
    trim_shiftb = [trim_shiftb(:,1),fliplr(trim_shiftb(:,2:end))].* maskb;
    trim_coeffb = [trim_coeffb(:,1),fliplr(trim_coeffb(:,2:end))].* maskb;
    %trim_shift = trim_shift .* maskb;
    %trimfold = max(cumsum(maskb,2),[],2);
    %figure;
    if plots == 1;  figure(91); imagesc(trim_shift); caxis([-5 5]); end;
    if plots == 1;  figure(90); imagesc(trim_shiftb); caxis([-5 5]); end;
    %imagesc(trim_shiftb);
    
    %trim_sum = sum(abs(trim_shift(:,startvol:endvol)),2) ./ max(cumsum(maskb(:,startvol:endvol),2),[],2);
    %     trim_shiftcj = trim_shift;
    %trim_datacj = trim_data;
    
    trim_shift = cumsum(trim_shift,2);
    trim_shift = conv2(filt_smo2,filt_smo,trim_shift,'same');
    trim_shift = trim_shift .* mask;
    %trim_shift(data==0) = 0;
    %         trim_sum = sum(trim_shift,2);
    %n = length(trim_shift);
    %trim_sum = sqrt((1/n)*sum((trim_shift.^2),2));
    
    
    %         % normalise the sum of the shifts for the fold in case we have
    %         % difference in fold
    %         norm_fold_b = max(cumsum(fmask(:,startvol:endvol),2),[],2);
    %         norm_fold = max(cumsum(fmask(:,startvol:midtrc),2),[],2);
    %         norm_fold_b(norm_fold_b == 0) = 1;
    %         norm_fold(norm_fold == 0) = 1;
    
    %         %tf   %stackin = sum(trim_data(:,startvol:endvol),2) ./ max(cumsum(maskb(:,startvol:endvol),2),[],2);
    %         stackin = sum(trim_data(:,startvol:midtrc),2) ./ norm_fold;
    %         stackinb = sum(trim_data(:,startvol:endvol),2) ./ norm_fold_b;
    %         %stackin = sum(trim_data(:,startvol:midtrc),2) ./ max(cumsum(fmask(:,startvol:midtrc),2),[],2);
    %         %stackinb = sum(trim_data(:,startvol:endvol),2) ./ max(cumsum(fmask(:,startvol:endvol),2),[],2);
    
    
    %if plots == 2; ilplotin(:,kk) = stackin; end
    ftrim_shift = interp1(time_sub,trim_shift,time(:,1),'linear',0);
    ftrim_shiftb = interp1(time_sub,trim_shiftb,time(:,1),'linear',0);
    
    trim_coeffb(isnan(trim_coeffb)) = 0;
    ftrim_coeffb = interp1(time_sub,trim_coeffb,time(:,1),'linear',0);
    maskpi = interp1(time_sub,mask,time(:,1),'linear',0);
    
    if strcmp(mode,'apply')
        for ii = 1:n_traces_gather
            out_traces(:,ii,kk) = interp1(time(:,ii),trim_data(:,ii),time(:,ii)-ftrim_shift(:,ii),'linear',0);
            trim_coeffb_out(:,ii) = interp1(time(:,ii),ftrim_coeffb(:,ii),time(:,ii)-ftrim_shift(:,ii),'linear',0);
        end
    else
        for ii = 1:n_traces_gather
            trim_coeffb_out(:,ii) = interp1(time(:,ii),ftrim_coeffb(:,ii),time(:,ii)-ftrim_shift(:,ii),'linear',0);
            out_traces(:,ii,kk) = interp1(time(:,ii),ftrim_shiftb(:,ii),time(:,ii)-ftrim_shift(:,ii),'linear',0);
        end
        
        out_traces(:,1:n_traces_gather,kk) = cumsum(out_traces(:,1:n_traces_gather,kk),2);
        out_traces(:,1:n_traces_gather,kk) = squeeze(out_traces(:,1:n_traces_gather,kk)) .* maskpi;
        

        
        
        % pick xcorrelation peaks - these are times to pick velocities at
        
        %             coeff_peaks_in = coeff_sum;
        %             coeff_peaks_out = zeros(size(coeff_peaks_in),1);
        
        
        %             while max(coeff_peaks_in)>0.95
        %
        %                 [~,peak_idx] = max(coeff_peaks_in);
        %                 coeff_peaks_out(peak_idx) = coeff_peaks_in(peak_idx);
        %                 start_mask=max(peak_idx-50,1);
        %                 end_mask=min(peak_idx+50,wb_n_samples);
        %
        %                 coeff_peaks_in(start_mask:end_mask) = 0;
        %
        %             end
        
    end
    %trim_coeffb_out(trim_coeffb_out == 0) = NaN;               % set zero to NaN for nanmean
    %coeff_sum(:,kk) = nanmean(trim_coeffb_out,2);
    %coeff_sum(isnan(coeff_sum)) = 0;
    
    trim_coeffb_out(isnan(trim_coeffb_out)) = 0;
    coeff_sum(:,kk) = mean(trim_coeffb_out,2);
    
    stack_traces = 0;
    stack_shifts=0;
    
    
    %         % %     % cj test version of shifts ========================
    %         %      zoom3 = 600;
    %         %      zoom4 = 1000;
    %         % %     figure(108); imagesc(trim_datacj(zoom3:zoom4,:)) ;  colormap(gray); caxis([-50 50]);
    %         % %     trim_datacjb = trim_datacj;
    %         % %     trim_datacj = fliplr(trim_datacj);
    %         % %     trim_shiftcj = fliplr(trim_shift);
    %         % %
    %         % % %      for ii = 2:n_traces_gather
    %         % % %          for jj = 1:ii
    %         % % %              trim_datacj(:,jj) = interp1(time(:,1),trim_datacj(:,jj),time(:,1)-trim_shiftcj(:,ii),'linear',0);
    %         % % %          end
    %         % % %      end
    %         % %
    %         % %      for ii = 2:n_traces_gather
    %         % %          trim_datacj(:,1:ii) = interp1(time(:,1),trim_datacj(:,1:ii),time(:,ii)-trim_shiftcj(:,ii),'linear',0);
    %         % %      end
    %         % %     trim_datacj = fliplr(trim_datacj);
    %         % %     %trim_shiftcj = fliplr(trim_shiftcj);
    %         %      figure(109); imagesc(trim_datacj(zoom3:zoom4,:)) ;  colormap(gray); caxis([-500000 500000]);
    %         %      figure(119); imagesc(trim_datacj) ;  colormap(gray); caxis([-500000 500000]);
    %         % %     figure(110); imagesc(trim_datacjb(zoom3:zoom4,:) - trim_datacj(zoom3:zoom4,:) ) ;  colormap(gray); caxis([-20 20]);
    %         % %     figure(112); imagesc(trim_datacjb(zoom3:zoom4,:) - trim_data(zoom3:zoom4,:) ) ;  colormap(gray); caxis([-20 20]);
    %         %      figure(121); imagesc(trim_data(zoom3:zoom4,:)) ;  colormap(gray); caxis([-500000 500000]);
    %         %      figure(111); imagesc(trim_data) ;  colormap(gray); caxis([-500000 500000]);
    %         %     % ===================================================
    %
    %
    %         if plots == 1;
    %             figure(102); imagesc(time_balence(trim_data)); colormap(gray); caxis([-4000 4000]);
    %             figure(103); imagesc(trim_data); colormap(gray);
    %         end;
    %
    %         % normalise the sum of the shifts for the fold in case we have
    %         % difference in fold
    % %         norm_fold_b = max(cumsum(fmask(:,startvol:endvol),2),[],2);
    % %         norm_fold = max(cumsum(fmask(:,startvol:midtrc),2),[],2);
    % %         norm_fold_b(norm_fold_b == 0) = 1;
    % %         norm_fold(norm_fold == 0) = 1;
    %
    %         %tf   %trim_sum = sum(abs(trim_shiftb_out(:,startvol:endvol)),2) ./ max(cumsum(maskb(:,startvol:endvol),2),[],2);
    %         trim_sum = sum(abs(trim_shiftb_out(:,startvol:endvol)),2) ./ norm_fold_b;
    %         %trim_sum = sum(abs(trim_shiftb_out(:,startvol:endvol)),2) ./ max(cumsum(fmask(:,startvol:endvol),2),[],2);
    %         %trim_sum = sum(abs(trim_shiftb_out(:,startvol:endvol)),2);
    %         if plots == 1;  figure(109); imagesc(trim_shiftb_out); caxis([-5 5]); end;
    %
    %         %tf   %stackout = sum(trim_data(:,startvol:endvol),2) ./ max(cumsum(maskb(:,startvol:endvol),2),[],2);
    %         stackout = sum(trim_data(:,startvol:midtrc),2) ./ norm_fold;
    %         %stackout = sum(trim_data(:,startvol:midtrc),2) ./ max(cumsum(fmask(:,startvol:midtrc),2),[],2);
    %         %stackout = sum(trim_data(:,startvol:endvol),2) ./ max(cumsum(fmask(:,startvol:endvol),2),[],2);
    %
    %         % =================================================================
    %         % now add a section to apply a residual trim static to the whole gather
    %         % based on the mismatch between the stacks before and after trim and
    %         % then repeat the final stack
    %
    %         st1 = double(time_balence_stk(stackin));
    %         sT1 = S3*spdiags(st1,0,wb_n_samples,wb_n_samples);
    %
    %         st2 = double(time_balence_stk(stackout));
    %         st2_shifts = ifft(bsxfun(@times,fft(st2),exp(1i*phase_shifts_orig)),'symmetric');
    %         %         for count = 1:totnumshifts
    %         %             sT2 = S3*spdiags(st2_shifts(:,count),0,wb_n_samples,wb_n_samples);
    %         %             det_coef5(:,count) =  ((sT1*st2_shifts(:,count)).*(sT2*st1))./((sT1*st1).*(sT2*st2_shifts(:,count)));
    %         %
    %         %         end
    %
    %
    %         countr = 1;
    %
    %         for count = 2:2:totnumshifts
    %             sT2 = S3*spdiags(st2_shifts(:,count),0,wb_n_samples,wb_n_samples);
    %             det_coef3(:,countr) =  ((sT1*st2_shifts(:,count)).*(sT2*st1))./((sT1*st1).*(sT2*st2_shifts(:,count)));
    %             countr = countr +1;
    %         end
    %         det_coef3(isnan(det_coef3)) = 0;
    %         det_coef4  = interp1(4:4:(totnumshifts*2),det_coef3',1:( ((totnumshifts-1)/2)*4 + 3  ),'spline');
    %         % add in a taper to zero out large shifts in the shallow before picking
    %         % the max
    %         det_coef4 = det_coef4'.*shift_maskstk;
    %         [~,stkidx]= max(det_coef4,[],2);
    %         stack_shift = (stkidx-(idxshift*2))*(shiftinc/2);
    %
    %         trim_data(:,1:n_traces_gather) = interp1(time(:,1),trim_data(:,1:n_traces_gather),time(:,1)+stack_shift(:,1),'linear',0);
    %
    %
    %         stackoutb = sum(trim_data(:,startvol:endvol),2) ./ norm_fold_b;
    %         stackout = sum(trim_data(:,startvol:midtrc),2) ./ norm_fold;
    %         %stackoutb = sum(trim_data(:,startvol:endvol),2) ./ max(cumsum(fmask(:,startvol:endvol),2),[],2);
    %         %stackout = sum(trim_data(:,startvol:midtrc),2) ./ max(cumsum(fmask(:,startvol:midtrc),2),[],2);
    %         % =================================================================
    %
    %         if plots == 2; ilplot(:,kk) = stackout; end
    %         %ilplotb(:,kk) = stackoutb;
    %
    %
    %
    %         results_out{2,2}(wb_idx(kk):wb_idx(kk)+wb_n_samples-1,:,kk) = ftrim_shift;
    %         traces(wb_idx(kk):wb_idx(kk)+wb_n_samples-1,:,kk) = trim_data;
    %
    %         results_out2{2,2}(wb_idx(kk):wb_idx(kk)+wb_n_samples-1,kk) = trim_sum;
    %         results_out2{3,2}(wb_idx(kk):wb_idx(kk)+wb_n_samples-1,kk) = stackin;
    %         results_out2{4,2}(wb_idx(kk):wb_idx(kk)+wb_n_samples-1,kk) = stackout;
    %         results_out2{5,2}(wb_idx(kk):wb_idx(kk)+wb_n_samples-1,kk) = stackoutb;
    %         results_out2{6,2}(wb_idx(kk):wb_idx(kk)+wb_n_samples-1,kk) = stackinb;
    %
    % %    end
    % end
    %
    % % if plots == 2;
    % %     zoom1 = 400;
    % %     zoom2 = 800;
    % %     %figure(102); imagesc((ilplotb(500:700,:) - ilplot(500:700,:))); colormap(gray);  caxis([-100 100])
    % %     %figure(101); imagesc(ilplotb(zoom1:zoom2,:)); colormap(gray);  caxis([-200 200])
    % %     figure(100); imagesc(ilplot(zoom1:zoom2,:)); colormap(gray);  caxis([-200 200])
    % %     figure(103); imagesc(ilplotin(zoom1:zoom2,:)); colormap(gray);  caxis([-200 200])
    % %     figure(102); imagesc((ilplotin(zoom1:zoom2,:) - ilplot(zoom1:zoom2,:))); colormap(gray);  caxis([-100 100])
    % %     %figure(104); imagesc((ilplotin(zoom1:zoom2,:) - ilplotb(zoom1:zoom2,:))); colormap(gray);  caxis([-100 100])
    % %     figure(105); imagesc(results_out2{2,2}(zoom1:zoom2,1:kk))
    % %     figure(107); imagesc(squeeze(results_out{2,2}(zoom1:zoom2,36,1:kk)))
    % %     figure(108); imagesc(reshape(results_out{3,2}(zoom1:zoom2,:,1:50:kk),(zoom2-zoom1+1),[]));  colormap(gray); caxis([-100 100]);
    % % end
    %
    % % need to reshape the 2 3d matricies into 2d ones as the segy write just wants samples * total traces
    %
    % results_out{2,2} = reshape(results_out{2,2},in_n_samples,[]);
    % traces = reshape(traces,in_n_samples,[]);
    %
    %
    %
    % % i_block = str2double(i_block);
    % % % write the pre-stack dataset
    % % node_segy_write(results_out,i_block, orig_sample_rate, output_dir)
    % % % write the stack
    % % node_segy_write(results_out2,i_block, orig_sample_rate, output_dir)
    
end