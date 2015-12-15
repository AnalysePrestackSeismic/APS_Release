function [] = trim_calculation_mute(job_meta_path,i_block,startvolr,endvolr,speedup,zsmooth,itmr,otmr,timbal,outputflatgathers,varargin)
%% ------------------ Disclaimer  ------------------
%
% BG Group plc or any of its respective subsidiaries, affiliates and
% associated companies (or by any of their respective officers, employees
% or agents) makes no representation or warranty, express or implied, in
% respect to the quality, accuracy or usefulness of this repository. The code
% is this repository is supplied with the explicit understanding and
% agreement of recipient that any action taken or expenditure made by
% recipient based on its examination, evaluation, interpretation or use is
% at its own risk and responsibility.
%
% No representation or warranty, express or implied, is or will be made in
% relation to the accuracy or completeness of the information in this
% repository and no responsibility or liability is or will be accepted by
% BG Group plc or any of its respective subsidiaries, affiliates and
% associated companies (or by any of their respective officers, employees
% or agents) in relation to it.
%% ------------------ License  ------------------
% GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
%% github
% https://github.com/AnalysePrestackSeismic/
%% ------------------ FUNCTION DEFINITION ---------------------------------
%% Function Description
%  This function works on migrated gathers to figure out the residual
%  statics through a sophisticated trim statics approach . It can flatten
%  the gather and also produce a gather of the prestack shifts required. It
%  can also make a measurement of average absolute value of shift required
%  across the gather to flatten the event which we can the MICA (Migrated
%  Image confidence Attrbute). The process has higher resolution than
%  conventional trim statics approac and more robust to the scale inssues
%  associated with trim statics


% ##################Data FLow Steps##################:
% I - Supply post mig input gathers
% II- Restrick offset with mutes: itmr and otmr
% III -Stack to produce [pre trim stack]*
% IV - Apply guide mute: startvolr,endvolr
% V -  Stack. This is [pre trim guide stack]*
% VI - Use Output(V) to generate the frequency dependednt parameters. Use this
% parameters on Output(II) to flatten the gather
% VII - Apply MICA algoithm to flatten Output(II)
% VIII - Measure MICA attribute as a stacked QC attribute. [trimsum]*
% IX - Apply guide mute and stack to produce a Guide Stack
% X - Use Guide stack to calculate bulkshifts to remove the jitter between
% Guide Stack and [pre trim guide stack].
% XI - Apply Output(X) to Output(VII) to produce [ Flattened Gathers]*
% XII - Apply guide mute to Output(XI) , stack to produce [post trim guide stack]*
% XIII- Stack Ouput(XI) to produce [post trim stack]*

% Output(Num) means Output after Step Num
% []* marked volumes are final outputs



% #####################Input:##############################
%   job_meta_path:      Path to job meta file
%   i_block:            The current Block Number
%   startvol:           The smallest angle trace to load to make the stack for stack alignment
%   endvol:             The largest angle trace to load to make the stack for stack alignment
%   speedup:            fast , or anything except fast
%   zsmooth :           Smoother ?? Not used?
%   itmr:               Inner Trace Mute for the trim zone, this is what is output in prestack gathers
%   otmr:               Outer Trace Mute for the trim zone, this is what is output in prestack gathers
%   timbal:             Flag to balence the data to avoid any loud parts dominating any solution. Use '1' for yes, '0' for no
%   outputflatgathers:  Flag to output the prestack flatterned data. Use '1' for yes, '0' for no
%   vargin:             Optional extra flags for zmin and zmax to limit the   trace length
% ######################Output:###########################
%   Void Output

% ###################### Writes to Disk:#######################
%   Everything gets written in directory : .../trimout/*.
%   The  volumes written out items are:
%       trimsum: Stack of the shifts (MICA attribute
%       pre trim guide stack: Stacks before final alignment of trimmed stack
%       post trim guide stack: Stacks after final alignment of trimmed stack
%       post trim stack: Stack of the trimmed gathers
%       pre trim stack: Stack of the untrimmed gathers

% Can also output depending on outputflatgathers
%       Flattened Gathers (after applying shifts)

% ########################Other Notes##########################
% QC
%   You can switch on various QC plots in the code assigning non zero values to plots variable
%   Check to code to understand what values to assign (1 or 2)

% For more details on MICA Algorithm please reference EAGE presentation or
% EAGE abstract . Ask Charles Jone, James, Selvage, Ratanadwip Ghosh or
% Dan Bright

% ============================================================================================================
%%
% sort out input parameters
samp_drop = 1; % set it to 2 if you want to drop every other sample. saves computation time
if strcmpi(speedup,'fast')  
    samp_drop = 2;
end
xcor_lowerlim = 0.85;
xcor_th  = 0.985; % for filtering trim shifts by an x-correlation threshold value
crap_data=0;
smoothing_scaler = 8;  % normally 10
startvolr = str2double(startvolr);
endvolr = str2double(endvolr);
otmr = str2double(otmr);
itmr = str2double(itmr);
timbal = str2double(timbal);
%shift = str2double(shift);
useselectemode = 0;
var_resample=0;
%outputflatgathers = 1; % 1 outputs the prestack flatterned data, 0 does not
outputflatgathers = str2double(outputflatgathers);

outputgathershifts = 0; % 1 outputs the prestack shifts, 0 does not
noofreports = 20;  % number of output reports in the log file

% pick up argements to truncate the z range
ztrunc = 0;
if ~isempty(varargin)
    % test to make sure zmin and zmax were set
    if length(varargin) == 2
        zmin = str2double(varargin{1});
        zmax = str2double(varargin{2});
        ztrunc = 1;
    else
        error('must specify zmin and zmax on command line if using z limits');
    end
end

%
if isempty(regexp(zsmooth,'il','once')) == 0
    useselectemode = 1;
    requiredinline =  str2double(regexprep(zsmooth,'il',''));
    %requiredinline =  str2double(strrep(tottracerun,'il',''));
    zsmooth = 20;
else
    zsmooth = str2double(zsmooth);
end

% zsmooth = str2double(zsmooth);
%zsmooth2 = shift*2;
% if zsmooth2 < 7
%     zsmooth2 = 7;
% end
% shiftinc = 0.5; now set in the code from the new sample rate
plots = 0;
padding = 50;
smxcorrlim = 0.87; % this is the threshold limit on the xcoor peaks, values less than this are set to 0
% Wavelet estimation parameters
ns_win = 128;
hns_win = ns_win/2;
ns_overlap = 96;
totmaxshift = 5;
avg_vel = 2500;
% ============================================================================================================
% load the job meta file containing more parameters
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
% old version but seismic structure is not used now
%[seismic, results_out{2,2}, results_out{1,2}{1,1}, offset_read] = node_segy_read(job_meta_path,'1',i_block);

% Read the data for this block
if ztrunc == 1
    [~, results_out{2,2}, results_out{1,2}{1,1}, offset_read] = node_segy_read(job_meta_path,'1',i_block,zmin,zmax);
else
    [~, results_out{2,2}, results_out{1,2}{1,1}, offset_read] = node_segy_read(job_meta_path,'1',i_block);
end
% should we update this for angle stacks?

% offset = unique(offset_read);
% results_out{2,2} = reshape(results_out{2,2},size(results_out{2,2},1),length(offset),[]);
%tabulate(offset_read);

offset_fold = size(results_out{2,2},2)/size(unique(results_out{1,2}{1,1},'rows'),1);

tmpoffset = tabulate(offset_read);
offset = tmpoffset(tmpoffset(:,2) ~= 0,1);
if offset_fold ~= size(offset,1)
    tmp2offset = offset([1;diff(offset,1)] > 1);
    offinc = mode(diff(offset([1;diff(offset,1)] > 1),1));
    if size(tmp2offset,1) < offset_fold
        offset = [tmp2offset;(((1:(offset_fold - size(tmp2offset,1))).*offinc)+offset(end))'];
    else
        offset = tmp2offset(1:offset_fold);
    end
    
end

results_out{2,2} = reshape(results_out{2,2},size(results_out{2,2},1),offset_fold,[]);

% put itmr and otmr into sequential trace numbers in the gather
itm = find(offset <= itmr,1,'last');
otm = find(offset <= otmr,1,'last');

% put startvolr and endvolr into sequential trace numbers in the gather and
% check if entries are valid
if (offset(1) <= startvolr)
    startvol = find(offset <= startvolr,1,'last');
    if (offset(end) > startvol)
        endvol = find(offset <= endvolr,1,'last');
    else
        fprintf('\n Provided endvolr invalid. This should be >= %d\n',startvol);
        error('Error: Wrong Input parameters');
    end
    
else
    fprintf('\n Provided startvolr invalid. This should be >= %d\n',offset(1));
    error('Error: Wrong Input parameters');
end


% check to make sure it read something if not exit
if isempty(results_out{2,2}) == 1 &&  isempty(results_out{1,2}{1,1}) == 1 && isempty(offset_read) == 1
    return
end

% drop data as required
% order is samples, angles, skeys
% results_out{2,2} = results_out{2,2}(1:maxzout,:,:);
% job_meta.n_samples{1} = maxzout;

%results_out{1,2}{1,1}

[n_samples,n_traces_gather,n_traces] = size(results_out{2,2});
in_n_samples = n_samples;

% apply any itm and otm to the gathers
sttrc = 2;
edtrc = n_traces_gather;
if itm > 0
    results_out{2,2}(:,1:itm,:) = 0;
    edtrc = (n_traces_gather - itm);
end
if otm > 0
    results_out{2,2}(:,otm:end,:) = 0;
    sttrc = (n_traces_gather - otm) + 2;
end
midtrc = startvol + floor((endvol - startvol)/2);
% ============================================================================================================
% Make results meta information

% for the pre-stack dataset
results_out{1,1} = 'Meta data for output files';
%results_out{resultno,2}{1,1} = ilxl_read;
results_out{1,2}{2,1} = offset_read';
ebcstrtowrite = sprintf('%-3200.3200s',[results_out{1,1} '  ' ebdichdr '  ' tmpebc]);
results_out{1,1} = ebcstrtowrite;
results_out{1,3} = 'is_gather'; % 1 is yes, 0 is no


% for the stack output
results_out2{1,1} = 'Meta data for output files';
results_out2{1,2}{1,1} = results_out{1,2}{1,1}(1:n_traces_gather:end,:);
results_out2{1,2}{2,1} = int32(zeros(n_traces,1));
%ebcstrtowrite = sprintf('%-3200.3200s',[results_out{resultno,1} '  ' ebdichdr '  ' tmpebc]);
results_out2{1,1} = ebcstrtowrite;
results_out2{1,3} = 'is_gather'; % 1 is yes, 0 is no

output_dir = [job_meta.output_dir,'trimout'];
% check to see if the directory exists
if exist(output_dir,'dir') == 0
    mkdir(output_dir);
end


% prestack output
filename = 'gaths';
results_out{2,1} = strcat(filename,'_trim_data_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
%results_out{2,2} = zeros(n_samples,n_traces_gather,n_traces,'single');
%commented out as written to above
results_out{2,3} = 1;

if outputgathershifts == 1;
    results_out{3,1} = strcat(filename,'_trim_shifts_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
    results_out{3,3} = 1;
    results_out{3,2} = zeros(in_n_samples,n_traces_gather,n_traces,'single');
end



% stack result
filename2 = 'stack';
results_out2{2,1} = strcat(filename2,'_mica_sum_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
results_out2{2,2} = zeros(in_n_samples,n_traces,'single');
results_out2{3,1} = strcat(filename2,'_mica_vals_good_data_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
results_out2{3,2} = zeros(in_n_samples,n_traces,'single');
%results_out2{4,1} = strcat(filename2,'_posttrim_no_resid_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
%results_out2{4,1} = strcat(filename2,'_posttrim_guide_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
results_out2{4,1} = strcat(filename2,'_mica_xcor_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
results_out2{4,2} = zeros(in_n_samples,n_traces,'single');
results_out2{5,1} = strcat(filename2,'_posttrim_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
results_out2{5,2} = zeros(in_n_samples,n_traces,'single');
results_out2{6,1} = strcat(filename2,'_pretrim_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
results_out2{6,2} = zeros(in_n_samples,n_traces,'single');
results_out2{2,3} = 0;
results_out2{3,3} = 0;
results_out2{4,3} = 0;
results_out2{5,3} = 0;
results_out2{6,3} = 0;
% ========================================================================
% flattern the data to the water bottom and make a stack of the data for
% freq analysis

% stack to make something to pick wb on
stackin_freq = sum(results_out{2,2}(:,startvol:midtrc,:),2) ./ ((midtrc - startvol) +1) ;
% pick wb

%% ==========================================================================================================================
% Pick water bottom or use a pre picked water bottom horizon
if isfield(job_meta, 'wb_path')
    wb_idx_in = dlmread(job_meta.wb_path);
    % col 1 inline
    % col 2 xline
    % col 3 twt
    if job_meta.is_gather == 1
        [~,locations] = ismember(results_out{1,2}{1,1}(1:length(offset):end,:),wb_idx_in(:,1:2),'rows');
    else
        error('should be gathers');
        %[~,locations] = ismember(ilxl_read{1}(1:end,:),wb_idx_in(:,1:2),'rows');
    end
    %wb_idx = zeros(size(traces{vol_index_wb},2),1);
    zero_loc = locations ~= 0;
    % make a zeros for each trace
    xi = (1:n_traces)';
    x = xi(zero_loc)';
    wb_idx = interp1(x,wb_idx_in(locations(zero_loc),3),xi);
    clear wb_idx_in
    
    wb_idx = (wb_idx./(job_meta.s_rate/1000))';
    
    wb_idx = round(wb_idx-padding);
    wb_idx(isnan(wb_idx)) = 1;
    wb_idx(wb_idx < 1) = 1;
    %win_sub = bsxfun(@plus,wb_idx,(0:job_meta.n_samples{1}-max(wb_idx))');
    
    %win_ind = bsxfun(@plus,win_sub,(0:job_meta.n_samples{1}:job_meta.n_samples{1}*(size(traces{vol_index_wb},2)-1)));
else
    %[wb_idx] = water_bottom_picker(stackin_freq,padding);
    [wb_idx] = water_bottom_picker(squeeze(stackin_freq),padding);
    wb_idx(wb_idx < 0) = 1;
end
% write out a wbpick
%ilxltoprint = results_out{1,2}{1,1}(1:length(offset):end,:);
% only write out the values that are not 1, so are picked
%dlmwrite(strcat(output_dir,'wbpick_',i_block,'.xyz'),[ilxltoprint((wb_idx ~= 1)',:),(int32(wb_idx(wb_idx ~= 1)+padding)'.*job_meta.s_rate)/1000],'delimiter', ' ', 'precision', '%-6d','newline', 'unix');
%
win_sub = bsxfun(@plus,wb_idx,(0:n_samples-max(wb_idx))');
win_ind = bsxfun(@plus,win_sub,(0:n_samples:n_samples*(n_traces-1)));
clear win_sub;
% make anew stack after flattening, make sure to blank it before
stackin_freq2 = stackin_freq(win_ind);
%stackin_freq2 = time_balence_stk(stackin_freq2);
[n_samples,~] = size(stackin_freq2);
clear win_ind stackin_freq;
%
%% ========================================================================
% find the dominant frequency of the data, take a subset of the data and
% average freqs down the trace
start_index = 1:ns_win-ns_overlap-1:n_samples-ns_win;
end_index = start_index+ns_win-1;
freq_axis = (1e6/job_meta.s_rate)/2*linspace(0,1,ns_win/2);
n_win = length(start_index);
peak_freq = zeros(n_win,2);
min_freq = zeros(n_win,2);
gathinc = ceil(n_traces/20);
% make taper to apply to signal before fft
taperlen = 16;
%taperst = linspace(0,1,taperlen)';
taperst = (sin(linspace((-pi/2),(pi/2),taperlen)')+1)/2;
taperend = 1 - taperst;
taperapply = [taperst;ones((ns_win-(taperlen*2)),1);taperend];
filt_smof =  ones(1,9)/9;
%
for ii = 1:n_win
    % Estimate wavelets and store meta information
    % make a tmp array with the fft values in it
    tmpfft = zeros(2,2);
    tmpfft = abs(fft(bsxfun(@times,stackin_freq2(start_index(ii):end_index(ii),:),taperapply)));
    avgfreq = sum(tmpfft,2,'double');
    avgfreqsmo = conv(avgfreq(1:hns_win),filt_smof,'same');
    %if plots == 1; figure(121); plot(avgfreqsmo); title('smoothed amplitude spectra');end;
    cja = (avgfreqsmo > max(avgfreqsmo)*0.75);
    [~,rhs] = max(cja(end:-1:1));
    [~,lowf_idx] = max(cja(1:1:end));
    highf_idx = (hns_win - rhs) + 1;
    %[~,idxf] = max(avgfreqsmo);
    %peak_freq(ii,:) = [(start_index(ii)+hns_win) ((peak_freq(ii,2)+freq_axis(highf_idx)))];
    peak_freq(ii,:) = [(start_index(ii)+hns_win) max(peak_freq(ii,2),freq_axis(highf_idx))];
    min_freq(ii,:) = [(start_index(ii)+hns_win) max(min_freq(ii,2),freq_axis(lowf_idx))];
    mid_freq(ii,:) = [(start_index(ii)+hns_win) (min_freq(ii,2)+((peak_freq(ii,2) -  min_freq(ii,2)  )*0.65))];
end
%
% =========================================================================
%see if the code can resample the data down to calculate shifts at
resamp_rt = round(1000 /(max(mid_freq(:,2))*2));
use_samp_drop = 0;
wb_n_samples = n_samples;
orig_sample_rate = sample_rate;
time = repmat((0:sample_rate:(n_samples-1)*sample_rate)',1,n_traces_gather);
if var_resample==1
    if resamp_rt > (sample_rate * 1.8)
        samp_drop = round(resamp_rt/(sample_rate*1.8) + 0.1);
        sample_rate = sample_rate * samp_drop;
        ufoff = 1000/(sample_rate*2);
        ufon = ufoff*0.8;
        n_samples = length(1:samp_drop:n_samples);
        peak_freq((peak_freq(:,2) > ufoff),2) = ufoff;
        totmaxshift = ceil(sample_rate/1.5); % was 2.5 and 3 before
        use_samp_drop = 1;
    end
else
    %     resamp_rt > (sample_rate * 1.8)
    %samp_drop = round(resamp_rt/(sample_rate*1.8) + 0.1); %set it to 2 if you want to drop every other sample. saves computation time
    sample_rate = sample_rate * samp_drop;
    ufoff = 1000/(sample_rate*2);
    ufon = ufoff*0.8;
    n_samples = length(1:samp_drop:n_samples);
    peak_freq((peak_freq(:,2) > ufoff),2) = ufoff;
    totmaxshift = ceil(sample_rate/1.5); % was 2.5 and 3 before
    use_samp_drop = 1;
end
time_sub = (0:sample_rate:(n_samples-1)*sample_rate);
%
% for jj = 1:gathinc:n_traces;
%     freq_test = results_out{2,2}(:,sttrc:midtrc,jj);
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
% ========================================================================
% calculate the shift limits and the smoothing for the trim, currently
% the wavelength at the peak freq at 2000m/s divided by 30
max_shift = [peak_freq(:,1) (avg_vel./(peak_freq(:,2)*42))];
% smooth the shifts
%filt_smos =  ones(1,3)/3;
%max_shift = conv(max_shift,filt_smos,'same');
% ========================================================================
% set the shift increment from the sample rate
shiftinc = sample_rate/30;  % this leads to floating point rounding errors
%shiftincloc = floor((sample_rate*10000)/30);


max_shiftstk = max_shift;
max_shiftstk(:,2) = max_shiftstk(:,2).*4; % was 3 before this allows the residual to do more

% clip the max shifts
max_shift((max_shift(:,2) > totmaxshift),2) = totmaxshift;
max_shiftstk((max_shiftstk(:,2) > totmaxshift),2) = (totmaxshift);

% set the maxshift to an increment of shiftinc
max_shift(:,2) = floor(max_shift(:,2)./shiftinc)*shiftinc;
max_shiftstk(:,2) = floor(max_shiftstk(:,2)./shiftinc)*shiftinc;

% calculate the length of the smoother, currently 10 times the shift value
freq_grid = [max_shift(:,1) floor(max_shift(:,2).*smoothing_scaler)];

%shift = max(max_shift(:,2));
shift = max(max_shiftstk(:,2));

%totmaxshift = shift;
% find the max shift to calculate
% if totmaxshift >= shift
%     totmaxshift = shift;
% else
%     shift = totmaxshift;
% end

% % clip the max shifts
% max_shift((max_shift(:,2) > totmaxshift),2) = totmaxshift;
% max_shiftstk((max_shiftstk(:,2) > totmaxshift),2) = (totmaxshift); % should multiple by 1.5 at some point, but need to update shift;

totnumshifts = length(-shift:shiftinc:shift);

%now make the time varying smoothing as a diagonal on a matrix S
%freq_grid2 = [30 11; 300 21; 900 31];

totalpts = max(freq_grid(:,2));
for mm = 1:size(freq_grid,1)
    freq_interp_grid(mm,:) = [zeros(1,(totalpts-freq_grid(mm,2))) , 0:(totalpts/freq_grid(mm,2)):totalpts ,  (totalpts-(totalpts/freq_grid(mm,2))):-(totalpts/freq_grid(mm,2)):0 , zeros(1,(totalpts-freq_grid(mm,2)))  ];
    
    %shift_interp_grid(mm,:) = [[zeros(1,round(shift/shiftinc-((max_shift(mm,2)/shiftinc)))) , ones(1,round(max_shift(mm,2)/shiftinc))],1,fliplr([zeros(1,round(shift/shiftinc-(max_shift(mm,2)/shiftinc))) , ones(1,round(max_shift(mm,2)/shiftinc))])];
    shift_interp_grid(mm,:) = [[zeros(1,round((round((shift/shiftinc)*10000) - round((max_shift(mm,2)/shiftinc)*10000))/10000)), ones(1,(round((max_shift(mm,2))/shiftinc)*10000)/10000)],1,fliplr([zeros(1,round((round((shift/shiftinc)*10000) - round((max_shift(mm,2)/shiftinc)*10000))/10000)),  ones(1,(round((max_shift(mm,2))/shiftinc)*10000)/10000)])];
    
    shift_interp_gridstk(mm,:) = [[zeros(1,round((round((shift/shiftinc)*10000) - round((max_shiftstk(mm,2)/shiftinc)*10000))/10000)), ones(1,(round((max_shiftstk(mm,2))/shiftinc)*10000)/10000)],1,fliplr([zeros(1,round((round((shift/shiftinc)*10000) - round((max_shiftstk(mm,2)/shiftinc)*10000))/10000)),  ones(1,(round((max_shiftstk(mm,2))/shiftinc)*10000)/10000)])];
    
    
    %shift_interp_grid(mm,:) = [[zeros(1,round(shift/shiftinc-((max_shift(mm,2)/shiftinc)))) , ones(1,round(max_shift(mm,2)/shiftinc))],1,fliplr([zeros(1,round(shift/shiftinc-(max_shift(mm,2)/shiftinc))) , ones(1,round(max_shift(mm,2)/shiftinc))])];
    %shift_interp_gridstk(mm,:) = [[zeros(1,round(shift/shiftinc-((max_shiftstk(mm,2)/shiftinc)))) , ones(1,round((max_shiftstk(mm,2))/shiftinc))],1,fliplr([zeros(1,round(shift/shiftinc-((max_shiftstk(mm,2))/shiftinc))) , ones(1,round((max_shiftstk(mm,2))/shiftinc))])];
    %shift_interp_grid(mm,:) = [[zeros(1,(shift/shiftinc-(max_shift(mm,2)/shiftinc))) , ones(1,max_shift(mm,2)/shiftinc)],1,fliplr([zeros(1,(shift/shiftinc-(max_shift(mm,2)/shiftinc))) , ones(1,max_shift(mm,2)/shiftinc)])];
    %shift_interp_gridstk(mm,:) = [[zeros(1,(shift/shiftinc-((max_shiftstk(mm,2)/shiftinc)))) , ones(1,(max_shiftstk(mm,2))/shiftinc)],1,fliplr([zeros(1,(shift/shiftinc-((max_shiftstk(mm,2))/shiftinc))) , ones(1,(max_shiftstk(mm,2))/shiftinc)])];
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

if plots == 1;
    figure(12); imagesc(diag_smo_interp); title('time varying windowing'); xlabel('number of samples'); ylabel('samples down trace')
    figure(13); imagesc(S); title('local windowing matrix for traces to go down the diagonal');
end;

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

% make a larger mask for small window xcoor
%shift_maskb = interp1(freq_zgrid,conv2(1,ones(1,(size(shift_interp_grid,2)/4))/(size(shift_interp_grid,2)/4),shift_interp_grid,'same'),start_interp:1:end_interp,'linear');
shift_maskb = interp1(freq_zgrid,conv2(1,ones(1,floor((size(floor(shift_interp_grid),2)/4))/floor((size((shift_interp_grid),2)/4))),shift_interp_grid,'same'),start_interp:1:end_interp,'linear');
shift_maskb(shift_maskb>0.1) = 1;
shift_maskb(shift_maskb<1) = 0;

shift_interp_gridstk = interp1(2:2:(totnumshifts*2),shift_interp_gridstk',1:((totnumshifts*2)+1),'linear');
shift_interp_gridstk(isnan(shift_interp_gridstk)) = 0;
shift_maskstk = interp1(freq_zgrid,shift_interp_gridstk',start_interp:1:wb_n_samples,'linear');
shift_mask(shift_mask<1) = 0;
shift_maskstk(shift_maskstk<1) = 0;
% ============================================================================================================
% loop round the traces 3d array one gather at a time and do calculations
f_max = (1/sample_rate)*1000;
phase_shifts = bsxfun(@times,(-shift:shiftinc:shift),(1/1000).*2.*pi.*repmat((0:f_max/(n_samples-1):f_max)',1,1+2*round((shift/shiftinc))));
if plots == 1; figure(14); imagesc(phase_shifts); title('phase shifts to apply to the trace'); colorbar('location','eastoutside'); end;
f_maxorig = (1/orig_sample_rate)*1000;
phase_shifts_orig = bsxfun(@times,(-shift:shiftinc:shift),(1/1000).*2.*pi.*repmat((0:f_maxorig/(wb_n_samples-1):f_maxorig)',1,1+2*round((shift/shiftinc))));

totnumshifts = length(-shift:shiftinc:shift);
%S = (1/zsmooth)*spdiags(repmat([(1:1:zsmooth),(zsmooth-1:-1:1)],n_samples,1),[(-zsmooth+1:1:0),(1:1:zsmooth-1)],n_samples,n_samples);
%S2 = (1/zsmooth2)*spdiags(repmat([(1:1:zsmooth2),(zsmooth2-1:-1:1)],n_samples,1),[(-zsmooth2+1:1:0),(1:1:zsmooth2-1)],n_samples,n_samples);
det_coef = zeros(n_samples,totnumshifts);
det_coefb = zeros(n_samples,totnumshifts);

%filttraces = [1 2 3 3 3 3 3 3 3 2 1]/27;
filttraces = [1 2 3 3 2 1]/12;
filt_smo =  [1 2 1]/4;
filt_smo2 =  ones(1,5)/5;
filt_smo3 =  [1 2 3 2 1]/9;
idxshift = round((shift/shiftinc))+1;

% initialise zeros in output
trim_shift = zeros(n_samples,n_traces_gather,'single');
trim_shiftb = zeros(n_samples,n_traces_gather,'single');
trim_xcorvb = zeros(n_samples,n_traces_gather,'single');
trim_shiftb_out = zeros(wb_n_samples,n_traces_gather,'single');
trim_xcorvb_out = zeros(wb_n_samples,n_traces_gather,'single');
wt_taperapply = [taperst;ones((wb_n_samples-(taperlen*2)),1);taperend];
det_coef3 = zeros(wb_n_samples,length(2:2:totnumshifts) );
norm_fold_b = zeros(wb_n_samples,1,'single');
norm_fold = zeros(wb_n_samples,1,'single');
%
close all;
reportinc = round(n_traces/noofreports);
%
for kk = 1:n_traces;
    
    if useselectemode == 0
        requiredinline = results_out{1,2}{1,1}(kk,1);
    end
    if results_out{1,2}{1,1}(kk,1) == requiredinline;
        
        
        % get the input gather and apply shift to flattern to the water bottom
        trim_data = results_out{2,2}(wb_idx(kk):(wb_idx(kk)+wb_n_samples-1),:,kk);
        
        %-----------------MUTES------------------------
        
        % apply an automatic mute to remove low amp noise
        fmask = low_amp_mute(trim_data);
        if isfield(job_meta,'mute') && kk==1
            mute_s=zeros(1,size(job_meta.mute,1));
            for m1=1:size(job_meta.mute,1)
               mute_s(m1)= find(offset <=job_meta.mute(m1,2) ,1,'last'); 
            end
            
            mute_func=floor(interp1(job_meta.mute(:,1),mute_s,1:size(trim_data,1)));
            mute_mask=zeros(size(trim_data ));
            for m2 = 1:size(trim_data,1)
                mute_mask(m2,1:mute_func(m2))=1;
            end
        trim_data = trim_data .* mute_mask;      
        end
        
        clear m1 m2;               
            
          
        trim_data = trim_data .* fmask;
        
        %--------------------------------------------------------
        
        if plots == 1; figure(106); imagesc(time_balence(trim_data)); colormap(gray); caxis([-4000 4000]); title('Input gather after scaling'); xlabel('offset/angle trace'); ylabel('samples down trace'); end;
        
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
        
        trim_data_filt_bk=trim_data_filt;
        if crap_data==1
              trim_data= medfilt3nt(trim_data,[0 5],0);
        end
        
        
        % balence the data to avoid any loud parts dominating any solution
        if max(max(trim_data_filt))==0 || isnan(max(max(trim_data_filt)))
            fprintf('Empty or Null traces present in block trace %d\n',kk);
            continue;
        end
        
        if timbal == 1
            trim_data_filtin = trim_data_filt;
            [trim_data_filt scaleused] = time_balence(trim_data_filtin);
            
            trim_data_filtin = trim_data_filt_bk;
            [trim_data_filt_bk ~] = time_balence(trim_data_filtin);
            clear trim_data_filtin;
        end
        
        %trim_data = trim_data_filt;
        %trim_data_filt = bandpass_filter(d,dt,f1,f2,f3,f4);
        
        if plots == 1; figure(101); imagesc(trim_data_filt); colormap(gray); caxis([-4000 4000]); title('Input gather after downsampling and scaling'); end; % Plot the input gather after downsampling
        tmp_plotin = trim_data_filt;
        % apply a small smoother to reduce noise
        if itm > 0
            mask(:,1:itm) = 0;
            trim_data_filt(:,1:itm) = repmat(trim_data_filt(:,(itm+1)),1,itm);
        end
        if otm > 0
            mask(:,otm:end) = 0;
            trim_data_filt(:,otm:n_traces_gather) = repmat(trim_data_filt(:,(otm-1)),1,(n_traces_gather-otm)+1);
        end
        
        % apply a small smoother to reduce noise
        %trim_data_filt = medfilt3nt(trim_data_filt,[0 9],0);
        trim_data_filt = conv2(filt_smo,filttraces,trim_data_filt,'same');
        
        
        
        %apply mask to the filtered data to remove duplicated data
        trim_data_filt = trim_data_filt .* mask;
        
        % flip the data to work from outside in
        trim_data_filt = fliplr(trim_data_filt);
        
        for ii = sttrc:edtrc
            t1 = double(trim_data_filt(:,ii-1));
            T1 = S*spdiags(t1,0,n_samples,n_samples);
            %             if plots == 1;
            %                 figure(15);  imagesc(T1); colormap(gray); title('Trace local windows');
            %                 figure(24); subplot(2,1,1); plot(t1); subplot(2,1,2); plot(diag(T1));title('compare input trace with trace down diagonal');
            %                 figure(25); subplot(2,1,1); plot(T1(200,:)); subplot(2,1,2); plot(T1(300,:)); title('local windows for xcorr');
            %             end;
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
            %det_coefb = det_coefb;
            %det_coefb = det_coefb.*shift_maskb;
            % cj edit took out the mask as it was making steps in the
            
            % results
            %            if plots == 1; figure(205); imagesc(det_coef); caxis([0.8 1]); end;
            %            if plots == 1; figure(206); imagesc(det_coefb); caxis([0.8 1]); end;
            %[~,idx] = max(det_coef');
            %trim_shift(:,ii-1) = idx'-shift-1;
            % this picks the index of the peak of the xcorelation idx and
            % the value of the xcorelation
            %[~,idx]= max(det_coef,[],2);
            [xcorva,idx]= max(det_coef,[],2);
            [xcorvb,idxb]= max(det_coefb,[],2);
            % set the threshold of data quality from the mean of the
            % xcorrelations on the whole trace
            xcor_th = (mean(det_coef(det_coef > 0.6))*0.995);
            xcor_thb = (mean(det_coefb(det_coefb.*shift_maskb > 0.6))*0.98);
            
            %            if plots == 1; figure(205); imagesc((-(idxshift-1):1:(idxshift-1))*shiftinc,1:1:size(det_coef,1),det_coef); caxis([0.8 1]); hold all; plot((idx-idxshift)*shiftinc,1:1:size(det_coef,1),'bs','MarkerEdgeColor','b','MarkerFaceColor','b'); xlabel('Shift'); ylabel('Sample No.'); end;
            %            if plots == 1; figure(206); imagesc((-(idxshift-1):1:(idxshift-1))*shiftinc,1:1:size(det_coefb,1),det_coefb); caxis([0.8 1]); hold all; plot((idxb-idxshift)*shiftinc,1:1:size(det_coefb,1),'bs','MarkerEdgeColor','b','MarkerFaceColor','b'); xlabel('Shift'); ylabel('Sample No.'); end;
            %            if plots == 1; figure(205); imagesc((-(idxshift-1):1:(idxshift-1))*shiftinc,1:1:450,det_coef(1:end,:)); caxis([0.8 1]); hold all; plot((idx(1:450)-idxshift)*shiftinc,1:1:450,'bs','MarkerEdgeColor','b','MarkerFaceColor','b'); xlabel('Shift'); ylabel('Sample No.'); title('smooth xcorrelation map'); end;
            %            if plots == 1; figure(206); imagesc((-(idxshift-1):1:(idxshift-1))*shiftinc,1:1:450,det_coefb(1:end,:)); caxis([0.8 1]); hold all; plot((idxb(1:450)-idxshift)*shiftinc,1:1:450,'bs','MarkerEdgeColor','b','MarkerFaceColor','b'); xlabel('Shift'); ylabel('Sample No.'); title('small smooth xcorrelation map'); end;
            
            %if plots == 1; figure; plot(idx,length(det_coef):-1:1); end;
            %if plots == 1; figure; plot(idxb,length(det_coef):-1:1); end;
 
            
            % insert the filtering of trimshifts by xcor to prevent uncorrelated event to jump too much            
            t_shift=(idx-idxshift)*shiftinc;                    % coarser trim shifts
%             subplot(3,1,2); plot(t_shift);grid on;
            t_shift_b=(idxb-idxshift)*shiftinc;                 % finer trim shifts

            % if data quality bad filter trim shift values by smoothened  xcor
            if crap_data==1
                xcorva_smooth= conv(xcorva,ones(50,1)/50,'same'); % Smoothen the xcor trace
                xcorvb_smooth= conv(xcorvb,ones(50,1)/50,'same'); % Smoothen the xcor trace
                t_shift(xcorva_smooth<xcor_th)=0;
%                 subplot(3,1,3); plot(t_shift);grid on;
                t_shift_b(xcorvb_smooth<xcor_thb)=0;
            else
                % now select low xcorr value areas and use with 0 to make a mask to
                % eliminate the shifts from those locations
                t_shift(isnan(xcorva))=0;
                t_shift_b(isnan(xcorvb))=0;               
                t_shift(xcorva<xcor_th)=0;
                t_shift_b(xcorvb<xcor_thb)=0;
            end
            
            trim_shift(:,ii) = t_shift;
            trim_shiftb(:,ii) = t_shift_b;            
            %trim_shift(:,ii) = (idx-idxshift)*shiftinc;
            %trim_shiftb(:,ii) = (idxb-idxshift)*shiftinc;
            trim_xcorvb(:,ii) = xcorvb;
            %trim_shift(:,ii-1) = idx-shift-1;
            
        end
        
        maskb = [ones(n_samples,1,'single'),mask(:,1:(end -1))] .* [mask(:,2:end),ones(n_samples,1,'single')];
        trim_shift = [trim_shift(:,1),fliplr(trim_shift(:,2:end))].* maskb;
        trim_shiftb = [trim_shiftb(:,1),fliplr(trim_shiftb(:,2:end))].* maskb;
        trim_xcorvb = [trim_xcorvb(:,1),fliplr(trim_xcorvb(:,2:end))].* maskb;
        
%         % now select low xcorr value areas and use with 0 to make a mask to
%         % eliminate the shifts from those locations
%         trim_xcorvb(isnan(trim_xcorvb))  = 0;
%         trim_xcorvb(trim_xcorvb<smxcorrlim) = 0;
%         trim_xcorvbmask = trim_xcorvb;
%         trim_xcorvbmask(trim_xcorvbmask>0) = 1;
%         trim_shiftb = trim_xcorvbmask.*trim_shiftb;
        
        % now scale the shift from z units to representing the average
        % wavelength, so using the same freq vertical distributin
        
        %trim_shift = trim_shift .* maskb;
        %trimfold = max(cumsum(maskb,2),[],2);
        %figure;
        if plots == 1;  figure(91); imagesc(trim_shift); caxis([-3 3]); title('gather of smooth shifts'); end;% Automate the indexing of trim_shift ???????????
        %if plots == 1;  figure(91); imagesc(trim_shift(1:450,1:25)); caxis([-3 3]); end;
        if plots == 1;  figure(90); imagesc(trim_shiftb); caxis([-5 5]); title('gather of small smooth shifts');end;
        if plots == 1;  figure(190); imagesc(trim_xcorvb); caxis([0.85 1]); title('gather of small smooth xcor peak values');end;
        %imagesc(trim_shiftb);
        %threshold the peak xcoors to select low xcor values and make a
        %mask to zero out those values in the sum
        
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
        
        
        % normalise the sum of the shifts for the fold in case we have
        % difference in fold
        norm_fold_b = max(cumsum(fmask(:,startvol:endvol),2),[],2);
        norm_fold = max(cumsum(fmask(:,startvol:midtrc),2),[],2);
        norm_fold_b(norm_fold_b == 0) = 1;
        norm_fold(norm_fold == 0) = 1;
        
        %tf   %stackin = sum(trim_data(:,startvol:endvol),2) ./ max(cumsum(maskb(:,startvol:endvol),2),[],2);
        stackin = sum(trim_data(:,startvol:midtrc),2) ./ norm_fold;
        stackinb = sum(trim_data(:,startvol:endvol),2) ./ norm_fold_b;
        %stackin = sum(trim_data(:,startvol:midtrc),2) ./ max(cumsum(fmask(:,startvol:midtrc),2),[],2);
        %stackinb = sum(trim_data(:,startvol:endvol),2) ./ max(cumsum(fmask(:,startvol:endvol),2),[],2);
        
        
        %if plots == 2; ilplotin(:,kk) = stackin; end
        ftrim_shift = interp1(time_sub,trim_shift,time(:,1),'linear',0);
        ftrim_shiftb = interp1(time_sub,trim_shiftb,time(:,1),'linear',0);
        ftrim_xcorvb = interp1(time_sub,trim_xcorvb,time(:,1),'linear',0);
        
        for ii = 1:n_traces_gather
            trim_data(:,ii) = interp1(time(:,ii),trim_data(:,ii),time(:,ii)-ftrim_shift(:,ii),'linear',0);
            %trim_data(:,ii) = interp1q(time(:,ii),trim_data(:,ii),time(:,ii)-trim_shift(:,ii));
            
            %trim_datab(:,ii) = interp1(time(:,ii),trim_data(:,ii),time(:,ii)-trim_shift(:,ii),'spline',0);
            trim_shiftb_out(:,ii) = interp1(time(:,ii),ftrim_shiftb(:,ii),time(:,ii)-ftrim_shift(:,ii),'linear',0);
            trim_xcorvb_out(:,ii) = interp1(time(:,ii),ftrim_xcorvb(:,ii),time(:,ii)-ftrim_shift(:,ii),'linear',0);
        end
        
        % %     % cj test version of shifts ========================
        %      zoom3 = 600;
        %      zoom4 = 1000;
        % %     figure(108); imagesc(trim_datacj(zoom3:zoom4,:)) ;  colormap(gray); caxis([-50 50]);
        % %     trim_datacjb = trim_datacj;
        % %     trim_datacj = fliplr(trim_datacj);
        % %     trim_shiftcj = fliplr(trim_shift);
        % %
        % % %      for ii = 2:n_traces_gather
        % % %          for jj = 1:ii
        % % %              trim_datacj(:,jj) = interp1(time(:,1),trim_datacj(:,jj),time(:,1)-trim_shiftcj(:,ii),'linear',0);
        % % %          end
        % % %      end
        % %
        % %      for ii = 2:n_traces_gather
        % %          trim_datacj(:,1:ii) = interp1(time(:,1),trim_datacj(:,1:ii),time(:,ii)-trim_shiftcj(:,ii),'linear',0);
        % %      end
        % %     trim_datacj = fliplr(trim_datacj);
        % %     %trim_shiftcj = fliplr(trim_shiftcj);
        %      figure(109); imagesc(trim_datacj(zoom3:zoom4,:)) ;  colormap(gray); caxis([-500000 500000]);
        %      figure(119); imagesc(trim_datacj) ;  colormap(gray); caxis([-500000 500000]);
        % %     figure(110); imagesc(trim_datacjb(zoom3:zoom4,:) - trim_datacj(zoom3:zoom4,:) ) ;  colormap(gray); caxis([-20 20]);
        % %     figure(112); imagesc(trim_datacjb(zoom3:zoom4,:) - trim_data(zoom3:zoom4,:) ) ;  colormap(gray); caxis([-20 20]);
        %      figure(121); imagesc(trim_data(zoom3:zoom4,:)) ;  colormap(gray); caxis([-500000 500000]);
        %      figure(111); imagesc(trim_data) ;  colormap(gray); caxis([-500000 500000]);
        %     % ===================================================
        
        
        if plots == 1;
            %figure(102); imagesc(time_balence(trim_data)); colormap(gray); caxis([-4000 4000]);
            tmp_plotout = time_balence(trim_data);
            plotinb = [fliplr(tmp_plotin),tmp_plotout(1:samp_drop:end,:)];% FLip left to right
            
            %figure(1020); imagesc(plotinb(1:450,:)); colormap(gray); caxis([-4000 4000]);
            
            figure(103); imagesc(tmp_plotout); colormap(gray); caxis([-4000 4000]); title('gather after raw smooth trim');
        end;
        
        % normalise the sum of the shifts for the fold in case we have
        % difference in fold
        %         norm_fold_b = max(cumsum(fmask(:,startvol:endvol),2),[],2);
        %         norm_fold = max(cumsum(fmask(:,startvol:midtrc),2),[],2);
        %         norm_fold_b(norm_fold_b == 0) = 1;
        %         norm_fold(norm_fold == 0) = 1;
        
        %tf   %trim_sum = sum(abs(trim_shiftb_out(:,startvol:endvol)),2) ./ max(cumsum(maskb(:,startvol:endvol),2),[],2);
        trim_sum = sum(abs(trim_shiftb_out(:,startvol:endvol)),2) ./ norm_fold_b;
        xcorvb_sum = sum(abs(trim_xcorvb_out(:,startvol:endvol)),2) ./ norm_fold_b;
        
        % make the mica_threshold version using a threshold on xcor to set
        % mica values to zero where data is poor quality
        %mica_threshold = trim_sum.* xcorvb_sum(xcorvb_sum < xcor_lo =  trim_sum;
        mica_threshold =  trim_sum;
        mica_threshold(xcorvb_sum < xcor_lowerlim) = 0;
        
        %trim_sum = sum(abs(trim_shiftb_out(:,startvol:endvol)),2) ./ max(cumsum(fmask(:,startvol:endvol),2),[],2);
        %trim_sum = sum(abs(trim_shiftb_out(:,startvol:endvol)),2);
        
        %if plots == 1;  figure(109); imagesc(trim_shiftb_out); caxis([-5 5]); end;
        if plots == 1;
            figure(109); imagesc(trim_shiftb_out); caxis([-5 5]); title('gather of small smooth shifts');
            figure(1091); imagesc(trim_xcorvb_out); caxis([0.85 1]);title('gather of small smooth xcor peak values');
            figure(1021); imagesc(plotinb); colormap(gray); caxis([-4000 4000]);  hold all;
            %plot(((trim_sum/max(trim_sum)).*2+9.5),1:1:length(trim_sum),'bs','MarkerEdgeColor','b','MarkerFaceColor','b');
            plot(((trim_sum/max(trim_sum)).*2+(size(plotinb,2)/2)),1:1:length(trim_sum),'bs','MarkerEdgeColor','b','MarkerFaceColor','b');
        end
        
        %tf   %stackout = sum(trim_data(:,startvol:endvol),2) ./ max(cumsum(maskb(:,startvol:endvol),2),[],2);
        stackout = sum(trim_data(:,startvol:midtrc),2) ./ norm_fold;
        %stackout = sum(trim_data(:,startvol:midtrc),2) ./ max(cumsum(fmask(:,startvol:midtrc),2),[],2);
        %stackout = sum(trim_data(:,startvol:endvol),2) ./ max(cumsum(fmask(:,startvol:endvol),2),[],2);
        
        % =================================================================
        % now add a section to apply a residual trim static to the whole gather
        % based on the mismatch between the stacks before and after trim and
        % then repeat the final stack
        
        
        if max(max(stackin))==0 || isnan(max(max(stackin)))
            fprintf('Empty or Null traces present in stackin trace %d\n',kk);
            continue;
        end
        
        st1 = double(time_balence_stk(stackin));
        
        sT1 = S3*spdiags(st1,0,wb_n_samples,wb_n_samples);
        
        if max(max(stackout))==0 || isnan(max(max(stackout)))
            fprintf('Empty or Null traces present in stackout trace %d\n',kk);
            continue;
        end
        st2 = double(time_balence_stk(stackout));
        
        st2_shifts = ifft(bsxfun(@times,fft(st2),exp(1i*phase_shifts_orig)),'symmetric');
        %         for count = 1:totnumshifts
        %             sT2 = S3*spdiags(st2_shifts(:,count),0,wb_n_samples,wb_n_samples);
        %             det_coef5(:,count) =  ((sT1*st2_shifts(:,count)).*(sT2*st1))./((sT1*st1).*(sT2*st2_shifts(:,count)));
        %
        %         end
        
        
        countr = 1;
        
        for count = 2:2:totnumshifts
            sT2 = S3*spdiags(st2_shifts(:,count),0,wb_n_samples,wb_n_samples);
            det_coef3(:,countr) =  ((sT1*st2_shifts(:,count)).*(sT2*st1))./((sT1*st1).*(sT2*st2_shifts(:,count)));
            countr = countr +1;
        end
        det_coef3(isnan(det_coef3)) = 0;
        det_coef4  = interp1(4:4:(totnumshifts*2),det_coef3',1:( ((totnumshifts-1)/2)*4 + 3  ),'spline');
        % add in a taper to zero out large shifts in the shallow before picking
        % the max
        det_coef4 = det_coef4'.*shift_maskstk;
        [~,stkidx]= max(det_coef4,[],2);
        stack_shift = (stkidx-(idxshift*2))*(shiftinc/2);
        
        trim_data(:,1:n_traces_gather) = interp1(time(:,1),trim_data(:,1:n_traces_gather),time(:,1)+stack_shift(:,1),'linear',0);
        if plots == 1;
            figure(1101);
            imagesc(time_balence(trim_data)); title('final outfut gather'); colormap('gray');caxis([-4000 4000]);
        end
        
        stackoutb = sum(trim_data(:,startvol:endvol),2) ./ norm_fold_b;
        stackout = sum(trim_data(:,startvol:midtrc),2) ./ norm_fold;
        %stackoutb = sum(trim_data(:,startvol:endvol),2) ./ max(cumsum(fmask(:,startvol:endvol),2),[],2);
        %stackout = sum(trim_data(:,startvol:midtrc),2) ./ max(cumsum(fmask(:,startvol:midtrc),2),[],2);
        % =================================================================
        
        
        % now make the mica values a function of the wavelength
        %using the peaks and toughs from the final stack stackoutb
        [sample_p, ~, ~, ~] = peak_picker(stackoutb);                      % pick peaks on the trace
        sample_p = sample_p*sample_rate;                                   % convert samples to time
        mica_wavelen = 100./([sample_p(1); (sample_p(2:end)-sample_p(1:end-1))]); % time between peaks to represent wavelength 100 div as preperation for %
        mica_wavelen = conv(mica_wavelen,filt_smo3,'same');                  % apply a smoother to remove jumps in wavelength
        mica_wavelength = interp1(sample_p,mica_wavelen,time(:,1),'linear',0);   % interp to the whole sample range
        mica_threshold = mica_threshold.*mica_wavelength;                       % calculate % mica shift rather than absolute shift
        
        if plots == 2; ilplot(:,kk) = xcorvb_sum; end
        %ilplotb(:,kk) = stackoutb;
        if outputgathershifts == 1;
            results_out{3,2}(wb_idx(kk):wb_idx(kk)+wb_n_samples-1,:,kk) = ftrim_shift; % Gather of the MICA attribute (finer shifts)
        end
        
        %-----------------FINAL RESULTS---------------------------
        results_out{2,2}(wb_idx(kk):wb_idx(kk)+wb_n_samples-1,:,kk) = trim_data;    % Gather with trim shifts applied
        
        results_out2{2,2}(wb_idx(kk):wb_idx(kk)+wb_n_samples-1,kk) = trim_sum;      % trimsum: Stack of the shifts (MICA attribute
        %results_out2{3,2}(wb_idx(kk):wb_idx(kk)+wb_n_samples-1,kk) = stackin;       % pre trim guide: Stacks before final alignment of trimmed stack
        results_out2{3,2}(wb_idx(kk):wb_idx(kk)+wb_n_samples-1,kk) = mica_threshold;       % mica values thresholded on xcor values
        results_out2{4,2}(wb_idx(kk):wb_idx(kk)+wb_n_samples-1,kk) = xcorvb_sum;      % Stack of the peakxcor values
        %results_out2{4,2}(wb_idx(kk):wb_idx(kk)+wb_n_samples-1,kk) = stackout;      % post trim guide: Stacks after final alignment of trimmed stack
        results_out2{5,2}(wb_idx(kk):wb_idx(kk)+wb_n_samples-1,kk) = stackoutb;     % post trim: Stack of the final trimmed gathers
        results_out2{6,2}(wb_idx(kk):wb_idx(kk)+wb_n_samples-1,kk) = stackinb;      % pre trim: Stack of the untrimmed gathers
        
        
        % Give a status report, add this to a log file
        
        curtr = kk/reportinc;
        curtrresid = curtr - round(curtr);
        if kk == 1;
            fprintf('Completed 0 percent: trace %d of %d \n',kk,n_traces)
        elseif (curtrresid == 0)
            fprintf('Completed %5.2f percent: trace %d of %d \n',((100/noofreports)*curtr),kk,n_traces)
        end

    end
end

if plots == 2;
    zoom1 = 400;
    zoom2 = 800;
    %figure(102); imagesc((ilplotb(500:700,:) - ilplot(500:700,:))); colormap(gray);  caxis([-100 100])
    %figure(101); imagesc(ilplotb(zoom1:zoom2,:)); colormap(gray);  caxis([-200 200])
    %figure(100); imagesc(ilplot(zoom1:zoom2,:)); colormap(gray);  caxis([-200 200])
    %figure(103); imagesc(ilplotin(zoom1:zoom2,:)); colormap(gray);  caxis([-200 200])
    %figure(102); imagesc((ilplotin(zoom1:zoom2,:) - ilplot(zoom1:zoom2,:))); colormap(gray);  caxis([-100 100])
    %figure(104); imagesc((ilplotin(zoom1:zoom2,:) - ilplotb(zoom1:zoom2,:))); colormap(gray);  caxis([-100 100])
    %figure(105); imagesc(results_out2{2,2}(zoom1:zoom2,1:kk))
    
    figure(106); imagesc(results_out2{5,2}(zoom1:zoom2,1:kk)); colormap(gray);  caxis([-200 200])
    figure(100); imagesc(ilplot(zoom1:zoom2,1:kk))
    figure(105); imagesc(results_out2{2,2}(zoom1:zoom2,1:kk))
    
    if outputgathershifts == 1;
        figure(107); imagesc(squeeze(results_out{3,2}(zoom1:zoom2,36,1:kk)))
    end
    figure(108); imagesc(reshape(results_out{2,2}(zoom1:zoom2,:,1:50:kk),(zoom2-zoom1+1),[]));  colormap(gray); caxis([-100 100]);
end

% need to reshape the 2 3d matricies into 2d ones as the segy write just wants samples * total traces
results_out{2,2} = reshape(results_out{2,2},in_n_samples,[]);
if outputgathershifts == 1;
    results_out{3,2} = reshape(results_out{3,2},in_n_samples,[]);
end




i_block = str2double(i_block);
% write the pre-stack dataset
if outputflatgathers == 1;
    if ztrunc == 1
        node_segy_write(results_out,i_block, orig_sample_rate, output_dir, zmin)
        %node_segy_write(results_out,i_block, orig_sample_rate, output_dir)
    else
        node_segy_write(results_out,i_block, orig_sample_rate, output_dir)
    end
end
% write the stack
if ztrunc == 1
    node_segy_write(results_out2,i_block, orig_sample_rate, output_dir, zmin)
    %node_segy_write(results_out2,i_block, orig_sample_rate, output_dir)
else
    node_segy_write(results_out2,i_block, orig_sample_rate, output_dir)
end
end

