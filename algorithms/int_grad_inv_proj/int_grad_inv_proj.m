function [] = int_grad_inv_proj(job_meta_path,i_block,startvol,volinc,endvol,tottracerun,maxzout,wavevar)
%
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
% INT_GRAD_INV_PROJ: function to run Intercept Gradient Inversion using
% dynamic wavelet set.
%   Inputs:
%       job_mat_path = path of metadata .mat file.
%       i_block = current block to be processed.
%       startvol = number of the first angle trace/volume to read
%       volinc = angle trace/volume increment
%       endvol = number of the last angle trace/volume to read
%       tottracerun = number of traces to output, 0 = all, and il1234 would
%       only output inline 1234
%       maxzout = To be given in samples if this is initialized as 0 this will be
%       converted to the maximum number of samples in the volume
%       wavevar = Flag 1: for spatially varying wavelet, Flag 0: for
%       spatially stationary wavelet

%   Outputs:
%       digi_intercept SEGY files.
%       digi_gradient SEGY files.
%       digi_minimum_energy_eer_projection SEGY files.
%
% -------------------------------------------------------------------------
%% Parameters
% would normaly convert all parameters to double, but keep i_block as string as being passed to
% other modules; it does happen at the bottom of this program for output
%i_block = str2double(i_block);
%
% angle trace data ranges to use, vol is an angle trace either as a
% seperate angle volume or as a angle trace in an angle gather

startvol = str2double(startvol);    % number of the first angle trace/volume to read
volinc = str2double(volinc);        % angle trace/volume increment 
%endvol = job_meta.nvols;
endvol = str2double(endvol);        % number of the last angle trace/volume to read

% number of traces to run, put to zero to make it run all traces in the
% block, this is the default, this is also used to pass an inline (pkey
% number to use in testing has to be ilnnnn format
useselectemode = 0;

if isempty(regexp(tottracerun,'il','once')) == 0
    useselectemode = 1;
    requiredinline =  str2double(regexprep(tottracerun,'il',''));
    %requiredinline =  str2double(strrep(tottracerun,'il',''));
    tottracerun = 0;
else
    tottracerun = str2double(tottracerun);
end
% tottracerun = 500;
%maxzout = 8000;
maxzout = str2double(maxzout);

output_std = 0;     % do you want to calculate and output standard line fitting intercept and gradient and eer
plot_on = 0;        % do you want interactive plots to pop up when running (just for debug)

% do you want to have a background starting model, default 0 is all zeros,
% setting this to value 1 uses standard line fitting to make a background
% model; lsq seems to reach the same answer with any model , did try random
% numbers and no difference in result just more iterations required
background = 0;
needconf = 0;                   % do you want a confidence volume calculating and outputing, 0 = no , 1 = yes
noofreports = 20;               % how many % progress reports do you want in the running of the job, 5 = 20% increment for reporting
chi_model_type = 'empirical';   % what is the chi model that you want for the eer...
                                % ...empirical is 19 degrees with 1 degree increase per km or second

% Tikhonov regularisation weight, 1 = little smoothing, 10 = moderate,
% 100 = smooth, 1000 = very smooth
%wsmooth = 7000; uruguay gathers
%wsmooth = 10; tza angle gathers
%wsmooth = str2double(wsmooth);

eer_weight = 0.1;   % Weight for EER constraint between 0 and 1, value of 1 forces it much  closer to the intercept values
iter = 300;         % maximum number of iterations in the lsq solver - normally 300
tol = 1e-3;         % convergence tolerance in the lsq solver can move to tol = 5e-4 or 1e-4.....but larger value speeds it up as less iterations to run 
%tol = 8e-4;  
padding = 64;       % water bottom pick padding
extrapad = 0;     % extra padding to top of dataset, this many samples will get zeroed below the wb 
%extrapad = (floor(extrapad/32))*32;
use_spatial_wavelets = wavevar; % FLag 1 for using spatial wavelets

warning off all;    % to reduce printout in compilied version turned all warning off

relfreq = 6;        % frequency to taper down below in hz, it is a sine taper so gentle ramp off, so can be a bit higher than expected, ie 6 cuts in at 4.5Hz
lowfreqtaperstart = 6; % Hz to start tapering the wavelet down from to 0 hz
wsmo_scal = 5; % scaler to make the weighting in the tikonov regularisation , normally set to 5 * the std dev of the amplitudes of the input data as measured
               % in the wavelet estimation code  

% end of parameters
%%
% set initial variables
totalvol = length(startvol:volinc:endvol);              % Total number of volumes to load
job_meta = load(job_meta_path);                         % Load job meta information 
wsmooth = job_meta.stdev_smo(str2double(i_block))*wsmo_scal;    % Get stdev for this block from the mat file
topfreq = 500000/job_meta.s_rate;
top3dpt = topfreq*0.72;
%
if tottracerun == 0
    eer_weight_out = num2str((eer_weight*1000));
    %eer_weight_out = regexprep(eer_weight_out, '0.', '');
    tol_out = num2str((tol*10000));
    %tol_out = regexprep(tol_out, '0.', '');
    
    ebdichdr = ['digi parameters: wsmooth ',num2str(wsmooth),' eer_weight ',num2str(eer_weight_out),' tolr ',tol_out];
    if useselectemode == 1;
        testdiscpt = ['_w_tolrt1e3__eer_weight_',num2str(eer_weight_out),'_il',num2str(requiredinline),'_range_',num2str(startvol),'_',num2str(volinc),'_',num2str(endvol)];
    else
        testdiscpt = ['_w_range_',num2str(startvol),'_',num2str(volinc),'_',num2str(endvol)];
    end    
else
    eer_weight_out = num2str((eer_weight*1000));
    %eer_weight_out = regexprep(eer_weight_out, '0.', '');
    tol_out = num2str((tol*10000));
    %tol_out = regexprep(tol_out, '0.', '');
    testdiscpt = [date,'_w_tolrt1e3__eer_weight_',num2str(eer_weight_out),'_range_p',num2str(startvol),'_',num2str(volinc),'_',num2str(endvol)];
    ebdichdr = ['digi parameters: wsmooth ',num2str(wsmooth),' eer_weight ',num2str(eer_weight_out),' tolr ',tol_out];
end

% add the history of jobs run and this one to the curent ebcdic
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

fprintf('reading data for total volumes = %d\n',totalvol)       % Make ouput directories and create meta information

% set maxzout if set to 0 on the command line
if maxzout == 0
    maxzout = job_meta.n_samples{1};                            % maxzout set yo maximum number of samples
end
%%
%==========================================================================================================================
% read data to find data to pick a water bottom on
vol_index_wb = 1;
if job_meta.is_gather == 0
    pick_wb_ind = ceil(job_meta.nvols*0.6667);
    
    [~, traces{vol_index_wb}, ilxl_read{vol_index_wb}] = node_segy_read(job_meta_path,num2str(pick_wb_ind),i_block);
    % check to make sure it read something if not exit
    if size(traces{vol_index_wb},1) == 1 && size(traces{vol_index_wb},2) == 1
        return
    end
    traces{vol_index_wb} = traces{vol_index_wb}(1:maxzout,:);
    n_samples_fullz = job_meta.n_samples{1};
    job_meta.n_samples{1} = maxzout;    
else
    
    
    % find the water bottom on a few gathers and take the middle one t use
    % as an offset plane to pick the water bottom on
    
    % read all the data for this block
    % node_segy_read(job_meta_path,vol_index,i_block)
    [~, vol_traces, ilxl_read{vol_index_wb}, offset_read] = node_segy_read(job_meta_path,'1',i_block);    
    vol_traces = vol_traces(1:maxzout,:);       % truncate data to make z axis number
    job_meta.n_samples{1} = maxzout;            % Write max zout in job meta file    
    offset = unique(offset_read);               % find the total number of offsets
    %#############for lobster TRN limit offset to 40 and rerun for failed
    %blocks###################################################################
    
    if isempty(vol_traces) == 1 &&  isempty(ilxl_read) == 1 && isempty(offset_read) == 1
        return
    end
    %juststack = 1;
    %if juststack == 1
    traces{vol_index_wb} = zeros(size(vol_traces,1),size(vol_traces,2)/size(offset,2));
    for stki  = 1:size(offset,2)
        traces{vol_index_wb} = traces{vol_index_wb} + vol_traces(:,offset_read == offset(stki));
    end
    
    pick_wb_ind = job_meta.live_offset_avg;
        
%         % should reshape and sum instead of this terrible loop
%     %else
%         
%         %read the middle 5 gathers in the input data to find the middle tkey(angle/offset) value of the water bottom
%         %tracestest{vol_index_wb} = vol_traces(:,(floor(size(vol_traces,2)/2)-(length(offset)*2)):(floor(size(vol_traces,2)/2)+(length(offset)*2)));
%         tracestest{vol_index_wb} = vol_traces(:,(floor(size(vol_traces,2)/2)-(length(offset)*2.5)+2):(floor(size(vol_traces,2)/2)+(length(offset)*2.5)));
%         tmpoffread = offset_read((floor(size(vol_traces,2)/2)-(length(offset)*2.5)+2):(floor(size(vol_traces,2)/2)+(length(offset)*2.5)));
%         %pick the water bottom
%         %pick the water bottom
%         [wb_idxcj] = water_bottom_picker(tracestest{vol_index_wb}(:,:),0);
%         %filter the water bottom pick to make a difference in WB time for
%         %traces that are not picking the wb
%         filtw = [1 2 3 2 1]/9;
%         wb_idxcjfilt = conv(wb_idxcj,filtw,'same');
%         % now make a blank array the size of the wb index array
%         wb_idxcj2 = zeros(1,size(wb_idxcj,2));
%         %now find the difference between the values of the filtered and
%         %unfiltered water bottom indexes, if there is a difference then it is
%         %likely to not be the water bottom as it should be mostly flat on the
%         %gathers
%         wb_idx_diff_ck = abs((wb_idxcj./wb_idxcjfilt)-1);
%         wb_idxcj2(wb_idx_diff_ck < 0.009) =  wb_idxcj(wb_idx_diff_ck < 0.009);
%         
%         %now work out the index locations of the water bottom
%         wb_idx_index = 1:1:size(wb_idxcj,2);
%         %apply a logical index to the index array to give the index of where the wb is picked and less then 10 elsewhere
%         wb_idx_index(ismember(wb_idxcj2,floor((min(wb_idxcj2((wb_idxcj2 > 10)))*0.9)):1:ceil((min(wb_idxcj2((wb_idxcj2 > 10)))*1.1))));
%         % select the angles with the wb on
%         tmpwbangs = (tmpoffread(wb_idx_index(ismember(wb_idxcj2,floor((min(wb_idxcj2((wb_idxcj2 > 10)))*0.9)):1:ceil((min(wb_idxcj2((wb_idxcj2 > 10)))*1.1))))));
%         tmpangpickstd = ceil(std(double(tmpwbangs)));
%         tmpangpick = floor(mean(tmpwbangs));
%         %tmpangpick = floor(mean(tmpoffread(wb_idx_index(ismember(wb_idxcj2,floor((min(wb_idxcj2((wb_idxcj2 > 10)))*0.9)):1:ceil((min(wb_idxcj2((wb_idxcj2 > 10)))*1.1)))))));
%         pick_wb_ind = find(offset == tmpangpick);
%         
%         clear tracestest offset_read;
% 
%         
% %         [wb_idxcj] = water_bottom_picker(traces{vol_index_wb}(:,:),0);
% %         
% %         %pick the water bottom
% %         [wb_idxcj] = water_bottom_picker(tracestest{vol_index_wb}(:,:),0);
% %         %filter the water bottom pick to make a difference in WB time for
% %         %traces that are not picking the wb
% %         filtw = [1 2 3 2 1]/9;
% %         wb_idxcjfilt = conv(wb_idxcj,filtw,'same');
% %         % now make a blank array the size of the wb index array
% %         wb_idxcj2 = zeros(1,size(wb_idxcj,2));
% %         %now find the difference between the values of the filtered and
% %         %unfiltered water bottom indexes, if there is a difference then it is
% %         %likely to not be the water bottom as it should be mostly flat on the
% %         %gathers
% %         wb_idx_diff_ck = abs((wb_idxcj./wb_idxcjfilt)-1);
% %         wb_idxcj2(wb_idx_diff_ck < 0.005) =  wb_idxcj(wb_idx_diff_ck < 0.005);
% %         
% %         %now work out the index locations of the water bottom
% %         wb_idx_index = 1:1:size(wb_idxcj,2);
% %         %apply a logical index to the index array to give the index of where the wb is picked and less then 10 elsewhere
% %         wb_idx_index(wb_idxcj2 == min(wb_idxcj2((wb_idxcj2 > 10))));
% %         % calculate the gather index in each gather by removing the integer
% %         % number of gathers from the index number and putting back to all being
% %         % the same angle index in each gather ie 1-46,1-46,1-46 etc... and
% %         % findingf the average value of all the indexes
% %         pick_wb_ind = floor(mean(((wb_idx_index(wb_idxcj2 == min(wb_idxcj2((wb_idxcj2 > 10)))))/length(offset) - floor((wb_idx_index(wb_idxcj2 == min(wb_idxcj2((wb_idxcj2 > 10)))))/length(offset)))*length(offset)));
% %         %read the offset plane from the input data gathers.
% %         %traces{vol_index_wb} = vol_traces(:,offset_read == offset(pick_wb_ind));
        
        
    %end
end


% pkey_inc_mode = mode(job_meta.pkey_inc);
% skey_inc_mode = mode(job_meta.skey_inc);
% 
% n_iline = (ilxl_read{vol_index_wb}(:,1)-min(ilxl_read{vol_index_wb}(:,1)))/pkey_inc_mode+1;
% n_xline = (ilxl_read{vol_index_wb}(:,2)-min(ilxl_read{vol_index_wb}(:,2)))/skey_inc_mode+1;
% skeyn = (max(ilxl_read{vol_index_wb}(:,2))-min(ilxl_read{vol_index_wb}(:,2)))/skey_inc_mode+1;
% pkeyn = (max(ilxl_read{vol_index_wb}(:,1))-min(ilxl_read{vol_index_wb}(:,1)))/pkey_inc_mode+1;
% lin_ind = ((n_iline-1).*skeyn)+n_xline;

%==========================================================================================================================
% Pick water bottom or use a pre picked water bottom horizon
if isfield(job_meta, 'wb_path')
    %wb_idx = dlmread(job_meta.wb_path,'delimiter','\t');
    wb_idx_in = dlmread(job_meta.wb_path);
    wb_idx_in = sortrows(wb_idx_in,[1 2]);
    % col 1 inline
    % col 2 xline
    % col 3 twt
    [~,locations] = ismember(ilxl_read{1}(1:end,:),wb_idx_in(:,1:2),'rows');    
    %wb_idx = zeros(size(traces{vol_index_wb},2),1);
    zero_loc = locations ~= 0;
    %wb_idx
    zero_loc = wb_idx_in(locations(zero_loc),3); 
    xi = (1:size(traces{vol_index_wb},2))';
    x = xi(zero_loc)';
    wb_idx = interp1(x,wb_idx_in(locations(zero_loc),3),xi);
    clear wb_idx_in
    %min_il = min(ilxl_read{vol_index_wb}(:,1));
    %max_il = max(ilxl_read{vol_index_wb}(:,1));
    %min_xl = min(ilxl_read{vol_index_wb}(:,2));
    %max_xl = max(ilxl_read{vol_index_wb}(:,2));
    %wb_idx2 = wb_idx(wb_idx(:,1) >= min_il & wb_idx(:,1) <= max_il & wb_idx(:,2) >= min_xl & wb_idx(:,2) <= max_xl,:);
    %    wb_idx = wb_idx(wb_idx(:,1) >= min_il & wb_idx(:,1) <= (max_il+1) & wb_idx(:,2) >= min_xl & wb_idx(:,2) <= (max_xl+1),:);
    %   wb_idx(:,3) = wb_idx(:,3)./(job_meta.s_rate/1000);
    wb_idx = (wb_idx./(job_meta.s_rate/1000))';
    
    wb_idx = round(wb_idx-padding);
    win_sub = bsxfun(@plus,wb_idx,(0:job_meta.n_samples{vol_index_wb}-max(wb_idx))');
    
    win_ind = bsxfun(@plus,win_sub,(0:job_meta.n_samples{vol_index_wb}:...
        job_meta.n_samples{vol_index_wb}*(size(traces{vol_index_wb},2)-1)));
else
    [wb_idx] = water_bottom_picker(traces{vol_index_wb},padding);
    wb_idx(wb_idx < 0) = 1;
    win_sub = bsxfun(@plus,wb_idx,(0:job_meta.n_samples{vol_index_wb}-max(wb_idx))');
                                                                      
    win_ind = bsxfun(@plus,win_sub,(0:job_meta.n_samples{vol_index_wb}:...
    job_meta.n_samples{vol_index_wb}*(size(traces{vol_index_wb},2)-1)));
   
end
%%
%==========================================================================================================================
% read in all the data and decimate as required (allready read in for gathers)
if job_meta.is_gather == 1
    % reshape the gather array to the same 3d matrix as the angle volumes
    %vol_tracescjtmp = reshape(vol_traces,size(traces{vol_index_wb},1),length(offset),size(traces{vol_index_wb},2));
    %vol_traces = vol_tracescjtmp(:,(startvol:volinc:endvol),:);
    %clear vol_tracescjtmp;
    % this can just be vol_traces_cjtmp is not needed, just put in to check
    % debug stages
    vol_traces = reshape(vol_traces,size(traces{vol_index_wb},1),length(offset),size(traces{vol_index_wb},2));
    
    % try applying a running average to the data to smooth out variations
    % on the angle gathers
        
    filttraces = [1 1 1]/3;
    %----- Finding the mute and discarding the points outside it by a mask----------------
    for ii = 1:1:size(vol_traces,3)
        % mute the low amplitude areas
        mask = low_amp_mute(vol_traces(:,:,ii));
        % apply a small smoother to reduce noise
        vol_traces(:,:,ii) = conv2(1,filttraces,vol_traces(:,:,ii),'same');
        vol_traces(:,:,ii) = vol_traces(:,:,ii) .* mask;
        %imagesc(vol_traces(:,:,ii));
        %colormap(gray);
        %caxis([-500 500]);
    end
    
    %drop the angles that are not needed
    % note the 3 trace average above to allow 2:1 decimation with a bit of
    % anti aliasing, not perfect but ...
    vol_traces = vol_traces(:,(startvol:volinc:endvol),:);    
    
    input_angles = double(offset(startvol:volinc:endvol));
    % resize the ilxl data to drop to stack fold
    tmp_ilxlrd = ilxl_read{vol_index_wb}(1:length(offset):end,:);
    ilxl_read{vol_index_wb} = tmp_ilxlrd;
    clear tmp_ilxlrd;
else
    % Load block for remaining angle traces
    % initialise the data array
    vol_traces = zeros(totalvol,n_samples_fullz,size(traces{vol_index_wb},2));
    
    % if tottracerun == 0;
    %     vol_traces = zeros(totalvol,size(traces{vol_index_wb},1),size(traces{vol_index_wb},2));
    % else
    %     vol_traces = zeros(totalvol,size(traces{vol_index_wb},1),tottracerun);
    % end

    vol_count = 1;
    for i_vol = startvol:volinc:endvol
        fprintf('reading data no of vol = %d\n',i_vol)

        % Read traces
        %[~, traces{vol_count}, ~, ~] = node_segy_read(job_meta_path,num2str(i_vol),i_block);
        [~, vol_traces(vol_count,:,:), ~, ~] = node_segy_read(job_meta_path,num2str(i_vol),i_block);

        % Flatten traces to water bottom
        %traces{vol_count} = traces{vol_count}(win_ind);
        %vol_traces(vol_count,:,:) = vol_traces(vol_count,win_ind);
        % truncate the dataset if just running a test for a set number of
        % traces
        %     if tottracerun ~= 0;
        % %         if vol_count == 1;
        % %             vol_traces = zeros(totalvol,size(traces{vol_count},1),tottracerun);
        % %         end
        %         %###########################################################
        %         % cj test edit resize the dataset temp to test with
        %         vol_traces(vol_count,:,:) = traces{vol_count}(:,1:tottracerun);
        %         traces{vol_count} = traces{vol_count}(:,1:tottracerun);
        %     end
        
        % keep a list of the angle values read in, this might need to change
        % for angle gathers with on a single angle value; it needs to be
        % assigned to both fields in the job meta file
        input_angles(vol_count) = (job_meta.angle{i_vol}(2)+job_meta.angle{i_vol}(1))/2;
        
        if plot_on == 2
            figure(1)
            subplot(1,totalvol,vol_count); imagesc(vol_traces(vol_count,:,:));
        end
        
        vol_count = vol_count + 1;
    end

    % deal with zeros in the IGmatrix build rather than as a mask as
    % already zeros, this is just to cover low amp noise
    %----- Finding the mute and discarding the points outside it by a mask----------------
    for ii = 1:1:size(vol_traces,3)
        % mute the low amplitude areas
        mask = low_amp_mute(vol_traces(:,:,ii)');
        % apply a small smoother to reduce noise
        %vol_traces(:,:,ii) = conv2(1,filttraces,vol_traces(:,:,ii),'same');
        vol_traces(:,:,ii) = vol_traces(:,:,ii) .* mask';
        
        %imagesc(vol_traces(:,:,ii));
        %colormap(gray);
        %caxis([-500 500]);
    end    

    % truncate data to make z axis number
    vol_traces = vol_traces(:,1:maxzout,:);
    job_meta.n_samples{1} = maxzout;
end
clear traces win_sub wb_idx;
%%
%==========================================================================================================================
% now flatten to the water bottom offset plane at a time
%
padding = padding + extrapad;
%
if job_meta.is_gather == 1
    for i_vol = 1:totalvol
        tmp_vol = squeeze(vol_traces(:,i_vol,:));
        vol_traces(1:size(win_ind,1),i_vol,:) = tmp_vol(win_ind);
    end
    % resize the traces to drop the trailing blanks after the shift to the wb
    vol_traces = vol_traces(1:size(win_ind,1),:,:);
    % taper the top 25 samples to avoid wb reflection creating artifacts due to
    % high amplitude

    % need to apply a taper that does not go to zero, make the wavelet try
    % to match zero and does not work
    %taper = single([ zeros(floor(padding-5),1)',linspace(0,1,20)])';
    %vol_traces(1:length(taper),:,:) = bsxfun(@times,vol_traces(1:length(taper),:,:),taper);
    
else
    for i_vol = 1:totalvol
        tmp_vol = squeeze(vol_traces(i_vol,:,:));
        % flattern to water bottom
        vol_traces(i_vol,1:size(win_ind,1),:) = tmp_vol(win_ind);
    end
    % resize the traces to drop the trailing blanks after the shift to the wb
    vol_traces = vol_traces(:,1:size(win_ind,1),:);
    % taper the top 25 samples to avoid wb reflection creating artifacts due to
    % high amplitude
    
    % need to apply a taper that does not go to zero, make the wavelet try
    % to match zero and does not work
    %taper = single([ zeros(floor(padding-5),1)',linspace(0,1,20)]);
    %vol_traces(:,1:length(taper),:) = bsxfun(@times,vol_traces(:,1:length(taper),:),taper);
    
end


ns = size(win_ind,1);
ntraces = size(win_ind,2);
%
% inverse taper to apply to the results in the inversion
%taperrev = single([ zeros(floor(padding-5),1)',linspace(1,0,20)]);
taperrev = single([ zeros(floor(padding-5),1)',1./(linspace(0,1,20)),ones((ns-(padding+15)),1)']);
taperrev = [taperrev, taperrev];
taperrev(isnan(taperrev)) = 0;
taperrev(isinf(taperrev)) = 0;
taperrev(taperrev> 10) = 10;
taperrev = taperrev';
clear tmp_vol;
padding = padding - extrapad;

% Test unflatten
% Unflatten data
% [ns,ntraces] = size(traces{1});
% traces_unflat = [traces{1};zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
% for kk = 1:length(wb_idx)
%     traces_unflat(:,kk) = circshift(traces_unflat(:,kk),wb_idx(kk));
% end
%
%==========================================================================================================================
%%
%---------Wavelets to be used-------------------------------------

%Load wavelets set
if use_spatial_wavelets == '0'                                                      % Use single wavelet set
    wavelets = load(strcat(job_meta.wav_directory,'all_wavelets_time.mat'));        % This file gets made by wavelet_average.m 
else                                                                                % Use Spatially Varying Wavelet Sets
    algo = '2';                                                                     % FLag for which algorithm to use for spatial wavelets 
    wavelet_dir_path = job_meta.wav_directory;                                      % Wavelet Directory Path from job meta file
    wavelets = load_wavelet_spatial(job_meta_path,wavelet_dir_path,i_block,algo);   % Call function to smoothen wavelet spatially by various algorithms and load the wavelet for this block
end
%%
ns_wavelet = size(wavelets.all_wavelets_time{1},1)-1;                               % Number of samples in a wavelet?
hns_wavelet = floor(ns_wavelet/2);                                                  % Half of number of samples?

%=======================================
%make a low freq taper for the wavelet to avoid blowing up noise

taperlen = floor((lowfreqtaperstart/(500000/job_meta.s_rate))*hns_wavelet)+1;
taperst = (sin(linspace((-(pi*0.94)/2),((pi*0.8)/2),taperlen)')+1)/2;
%taperst(1:(end-2)) = taperst(1:(end-2)) ./ 10;
taperst= power(taperst,4);
taperend = flipud(taperst);
taperapply = [taperst;ones((ns_wavelet-(taperlen*2)),1);taperend];

%==========================================

vol_count = 1;
%Loop though the different angle volumes
for i_vol = startvol:volinc:endvol
    wavelet_z_grid = wavelets.all_wavelets_freq{i_vol}(1,:);
    wavelet_z_grid = wavelet_z_grid + padding;
    wavelet = wavelets.all_wavelets_freq{i_vol}(2:end,:);
    wavelet(isnan(wavelet)) = 0;
    % apply the low freq taper
    wavelet = bsxfun(@times,wavelet,taperapply);
    if plot_on == 1
        figure(2)
        subplot(1,totalvol,vol_count); imagesc(wavelet);
        figure(3)
        subplot(1,totalvol,vol_count); plot(wavelet(:,10));
    end
    start_interp = min(wavelet_z_grid);%-hns_wavelet;
    end_interp = max(wavelet_z_grid);%+hns_wavelet;
    wavelet_interp{vol_count} = interp1(wavelet_z_grid,wavelet',start_interp:1:end_interp,'linear');
    %wavelet_interp{vol_count} = interpolate_wavelets(wavelet_z_grid,wavelet,start_interp,end_interp);
    wavelet_interp{vol_count} = circshift(ifft(wavelet_interp{vol_count}','symmetric'),floor(job_meta.ns_win/2));
    %wavelet_interp{i_vol} = wavelet_interp{i_vol};
    
    if start_interp > 1;
        % Pad with zeros or first wavelet
        pad_size = start_interp-1;
        wavelet_interp{vol_count} = [repmat(wavelet_interp{vol_count}(:,1),1,pad_size),...
            wavelet_interp{vol_count}];
    end
    if end_interp < ns;
        % Pad with zeros or last wavelet
        pad_size = ns-end_interp;
        wavelet_interp{vol_count} = [wavelet_interp{vol_count},...
            repmat(wavelet_interp{vol_count}(:,end),1,pad_size)];
    end
    if floor(ns_wavelet/2) == ns_wavelet/2
        wavelet_interp{vol_count}(end+1,:) = wavelet_interp{vol_count}(end,:);
    end
    if plot_on == 1
        figure(4)
        subplot(1,totalvol,vol_count); imagesc(wavelet_interp{vol_count});
    end
    wavelet_interp{vol_count} = wavelet_interp{vol_count}(:,1:ns);
    vol_count = vol_count + 1;    
end
interp_wavelet_z_grid = 1:ns;
ns_wavelet = size(wavelet_interp{1},1);

%==========================================================================================================================
% normalise the wavelets
if job_meta.is_gather == 0
    [~, wbliveidx] =  min(abs((startvol:volinc:endvol)-pick_wb_ind));
    %wavelet_norm = wavelet_rms_norm(totalvol,wavelet_interp,interp_wavelet_z_grid,ceil(totalvol*0.6667));
    wavelet_norm = wavelet_rms_norm(totalvol,wavelet_interp,wbliveidx);
else
    % need to store the water bottom live offset in the job_meta
    %job_meta.livewb = 12;
    %[~, wbliveidx] = min(abs(input_angles-job_meta.livewb));
    [~, wbliveidx] =  min(abs((startvol:volinc:endvol)-pick_wb_ind));
    wavelet_norm = wavelet_rms_norm(totalvol,wavelet_interp,wbliveidx);
end

clear wavelet_interp wavelets wavelet_z_grid interp_wavelet_z_grid wavelet ;

for i_vol = 1:1:totalvol
    wavelet_norm{i_vol}(isnan(wavelet_norm{i_vol})) = 0;
end

if plot_on == 1
    figure(5)
    waveletck = cell2mat(wavelet_norm);
    imagesc(waveletck);
    imagesc(waveletck(:,1000:ns:end));
    plot(sum(abs(waveletck(:,1600:ns:end)),1));
end


%%
% start building the inversion operators
% Chi model
switch chi_model_type
    case 'empirical' 
        chi = (job_meta.s_rate/1e6)*(0:1:ns-1)'.*-2 + 19;

    case 'raw'
        
    case 'bootstrap'

    otherwise
        warning('Unexpected plot type. No plot created.');
end

% Build operator matrix
% Build blanking matrix used to ensure the convolution operator matrix is correct
IGblank = spdiags(ones(ns,2*hns_wavelet+1),(-hns_wavelet:hns_wavelet),ns,ns);
IGblank = repmat(IGblank,1+totalvol,2);

% Tikhonov regularisation matrix % check tik order
%##smooth = spdiags([-wsmooth*ones(2*ns,1) wsmooth*ones(2*ns,1)],[-1 0],2*ns,2*ns);
%smooth = spdiags([-wsmooth*ones(2*ns,1) wsmooth*ones(2*ns,1)],[-1 0],ns,ns);

smooth = spdiags([wsmooth*ones(2*ns,1) -2*wsmooth*ones(2*ns,1) wsmooth*ones(2*ns,1)],[-1 0 1],2*ns,2*ns);
%smooth = smooth.*0;
%cj edit
%smooth = spdiags([wsmooth*ones(2*ns,1) -2*wsmooth*ones(2*ns,1) wsmooth*ones(2*ns,1)],[-1 0 1],ns,ns);
%smooth2 = spdiags([wsmooth*10*ones(2*ns,1) -2*wsmooth*10*ones(2*ns,1) wsmooth*10*ones(2*ns,1)],[-1 0 1],ns,ns);
% smooth = [smooth, smooth];

% Extract the angle of the stack closest to the mean chi angle for this trace
% This is for the IG crossplot correlation constraint
theta = asind(sqrt(tand(mean(chi))));
[~, angle_index] = min(abs(input_angles-theta));

IGmatrix = build_operator(totalvol,input_angles,ns_wavelet,wavelet_norm,ns,hns_wavelet,angle_index,chi,eer_weight,IGblank,smooth);

clear wavelet_norm IGblank smooth ;


%% Inversion loop

% Set first and last traces to loop over
first_iter = 1;
if tottracerun ~= 0;
    ntraces = tottracerun;
end
last_iter = ntraces;

% Begin inversion loop
tic
ava = zeros(2*ns,((last_iter - first_iter)+1),'single');
%residcount = zeros(iter,(last_iter - first_iter)+1);
%cjmodel = rand(4100,1);
%
% make a taper to apply before any filterinf
tflen = 0.05;
%zeropad = 5;
up = ((sin(linspace((-(pi)/2),((pi)/2),round(tflen*ns))')+1)/2);
%down = flipud(up);
%outtaper = [up; ones(ns-(length(up)*2)-zeropad,1); down; zeros(zeropad,1); up; ones(ns-(length(up)*2)-zeropad,1); down; zeros(zeropad,1) ];
down = ((sin(linspace(((pi)/2),(-(pi)/2),round((tflen*2)*ns))')+1)/2);
outtaper = [up; ones(ns-(length(up)+length(down)),1); down; up; ones(ns-(length(up)+length(down)),1); down ];


% now loop round all the traces in the input and run the inversion we have
% added a parfor loop for testing but will not compile
for kk = first_iter:last_iter 
%parfor kk = first_iter:last_iter
%     for ii = 1:totalvol
%         data(:,ii) = traces{ii}(:,kk); % Read the angle stack data for the inversion of this trace
%     end

if useselectemode == 0
    requiredinline = ilxl_read{vol_index_wb}(kk,1); 
end

if ilxl_read{vol_index_wb}(kk,1) == requiredinline;

    if job_meta.is_gather == 1
        data = vol_traces(:,:,kk);
    else
        data = vol_traces(:,:,kk)';
    end
    data(isnan(data)) = 0;  % Set NaNs to zero
      
    fold = sum(data ~= 0,2); % Get angle fold  
    
    fmask = data == 0;
    %fmask = low_amp_mute(trim_data);

    % filter the low freq out of the dataset
        %[dataqc scaleused] = time_balence(data);
        %figure(57); imagesc(data); colormap(gray); caxis([-10000 10000]); 
        %amp_spec_plot(data,job_meta.s_rate)
    data = bandpass_filter(data,(job_meta.s_rate/1000000),0,relfreq,top3dpt,topfreq);
        %amp_spec_plot(data,job_meta.s_rate)
        %[dataqc scaleused] = time_balence(dataqc);
        %figure(58); imagesc(data); colormap(gray); caxis([-10000 10000]); 
    data = data .* (1.-fmask);
    %figure(59); imagesc(data); colormap(gray);
    data_tmp = data(:);  % Make temporary column vector from data
    
    %data_zeros = abs(data_tmp) < abs(mean(data_tmp)/range(data_tmp)); % Find zones where data is zero (due to mute angle mute functions)
    data_zeros = data_tmp == 0;
    data_zeros2 = data(:,wbliveidx) == 0;
    % or close to zero
    %data(reshape(data_zeros,[],totalvol)) = 0;
    %data_zeros = logical([data_zeros;zeros(3*ns,1)]); 
    %data_zeros = logical([data_zeros;data_zeros2;zeros(2*ns,1)]);
    data_zeros = logical([data_zeros;data_zeros2;data_zeros2;data_zeros2]); 
    
    if plot_on == 1;
        figure(6)
        subplot(1,2,1); spy(IGmatrix)
    end
    IGiter = IGmatrix;
    IGiter(data_zeros,:) = 0; % Set operator rows to zero if there are zeros in the data vector
    
    if plot_on == 1;
        figure(6)
        subplot(1,2,2); spy(IGiter)
    end
    
    % Model intercept = [near], model gradient = [far-near]/[far_angle - near_angle]
    if background == 1 || output_std == 1;
        model_tmp = zeros(2,ns);
        for ii = 1:ns
            model_op = [ones(totalvol,1),sind(input_angles').*sind(input_angles')];
            model_zeros = data(ii,:) == 0;
            model_op(model_zeros,:) = 0;
            model_tmp(:,ii) = model_op\data(ii,:)';
        end
        model_tmp = model_tmp';
        %[traces{1}(:,1) traces{2}(:,1) traces{3}(:,1) model_tmp]; 
        if output_std == 1;
            Imodel(:,kk) = model_tmp(:,1)/norm(model_tmp(:,1)); %data(:,1)/norm(data(:,1));
            Gmodel(:,kk) = model_tmp(:,2)/norm(model_tmp(:,1)); %-Imodel./tand(chi);
            model = [Imodel(:,kk);Gmodel(:,kk)];
            model(isnan(model)) = 0;
        else
            Imodel = model_tmp(:,1)/norm(model_tmp(:,1)); %data(:,1)/norm(data(:,1));
            Gmodel = model_tmp(:,2)/norm(model_tmp(:,1)); %-Imodel./tand(chi);
            %model = [Imodel;Gmodel];
            model(isnan(model)) = 0;
        end
    end
    % Set NaNs to zero
     
    % Make the data a column vector and add zeros on the end for the EER constraint and the Tikhonov regularisation
    data = double([data(:);zeros(3*ns,1,'single')]);
    %data = [data(:);zeros(2*ns,1)];
    % Do the inversion
    if background == 1;
        [ava_tmp,lsqflag,~,fiternotmp] = lsqr(IGiter,data,tol,iter,[],[],model);
        ava(:,kk) = single(ava_tmp);
        % try a different solving method?
        ava(:,kk) = bsxfun(@times,cjava,taperrev);
    else     
        %[ava(:,kk),~] = lsqr(IGiter,data,tol,iter,[],[],cjmodel);
        %[ava_tmp,lsqflag,relres,fiternotmp,residvec] = lsqr(IGiter,data,tol,iter,[],[]);
        [ava_tmp,lsqflag,~,fiternotmp] = lsqr(IGiter,data,tol,iter,[],[]);
        %residcount(1:length(residvec),kk) = residvec;
        %apply band pass filter to the output
        ava_zeros = ava_tmp == 0;
        % need to apply tapering to the data before the bandpass filter at
        % the start of each cube
        ava_tmp = ava_tmp.*outtaper;
        ava_tmp = bandpass_filter(ava_tmp,(job_meta.s_rate/1000000),(relfreq/3),(relfreq+1),top3dpt,topfreq);
        ava(:,kk) = ava_tmp .* (1.-ava_zeros);
        ava(1,kk) = single(fiternotmp);
        % left out the taper reverse to see if that just kept high
        % amplitudes down
        %ava(:,kk) = bsxfun(@times,ava_tmp,taperrev);
    end
    %[ava(:,kk),~] = lsqr(IGiter,data,1e-2,100,[],[],model);
    
    % mute out any hard zeros from the data
    %ava([fold;fold]==0,kk)=0;
    
%     % apply band pass filter to the output 
%     ava_zeros = ava(:,kk) == 0;
%     ava(:,kk) = bandpass_filter(ava(:,kk),(job_meta.s_rate/1000000),0,relfreq,top3dpt,topfreq);
%     ava(:,kk) = ava(:,kk) .* (1.-ava_zeros);
    
    % Estimate the R^2 confidence in the result. This is the variance ratio:
    % 1-[(sum(data-Gm)^2)/(sum(data-mean)^2)], where the inversion was
    % solving data = Gm.
    if needconf == 1;
        data = reshape(data(1:ns*totalvol,:),[],totalvol);
        digi_confidence(:,kk) = 1-(sum((data-reshape(IGiter(1:totalvol*ns,:)*ava(:,kk),[],totalvol)).^2,2)./sum(bsxfun(@minus,data,sum(data,2)./fold).^2,2));
    end
    % Give a status report, add this to a log file
    reportinc = round((last_iter-first_iter+1)/noofreports);
    curtr = (kk-first_iter+1)/reportinc;
    curtrresid = curtr - round(curtr); 
    if kk == first_iter;
        fprintf('Completed 0 percent: trace %d of %d lsqflag was: %d it took: %d iterations\n',kk-first_iter+1,last_iter-first_iter+1,lsqflag,fiternotmp)
    elseif (curtrresid == 0) 
        fprintf('Completed %5.2f percent: trace %d of %d lsqflag was: %d it took: %d iterations\n',((100/noofreports)*curtr),kk-first_iter+1,last_iter-first_iter+1,lsqflag,fiternotmp)
    end
end    
end
toc
%fprintf('Completed %8d traces: lsqflag was: %d last iteration no was: %d \n',(kk,lsqflag,fiternotmp);
clear IGiter IGmatrix;      % remove the preconditioning matrices

%%
% Save results
%digi_intercept = zeros(job_meta.n_samples{vol_index_wb},ntraces);
%digi_gradient = zeros(job_meta.n_samples{vol_index_wb},ntraces);

%digi_intercept = [ava(1:ns,:);zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
%digi_gradient = [ava(1+ns:end,:);zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
%digi_minimum_energy_eer_projection = [bsxfun(@times,ava(1:ns,:),cosd(chi))+bsxfun(@times,ava(1+ns:end,:),sind(chi));zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
%digi_minimum_energy_eer_projection = bsxfun(@times,ava(1:ns,:),cosd(chi))+bsxfun(@times,ava(1+ns:end,:),sind(chi));


if needconf == 1;
    digi_confidence(digi_confidence<0)=0;
    digi_confidence(digi_confidence>1)=1;
    digi_confidence(isnan(digi_confidence))=0;
    digi_confidence = [digi_confidence;zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
end
if output_std == 1;
    std_intercept = [Imodel;zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
    std_gradient = [Gmodel;zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
    %calculate the mimum energy based on chi angle provided
    std_minimum_energy_eer_projection = [bsxfun(@times,Imodel,cosd(chi))+bsxfun(@times,Gmodel,sind(chi));zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
end

% Old way to Unflatten data
%digi_intercept(win_ind(:,1:ntraces)) = ava(1:ns,:);
%digi_gradient(win_ind(:,1:ntraces)) = ava(1+ns:end,:);
%digi_minimum_energy_eer_projection = bsxfun(@times,digi_intercept,cosd(chi))+bsxfun(@times,digi_gradient,sind(chi));
%digi_confidence(win_ind(:,1:ntraces)) = digi_confidence(1:ns,:);
% for kk = 1:ntraces
% %cj edit    
% %for kk = 1:length(wb_idx)
%     digi_intercept(:,kk) = circshift(digi_intercept(:,kk),wb_idx(kk));
%     digi_gradient(:,kk) = circshift(digi_gradient(:,kk),wb_idx(kk));
%     digi_minimum_energy_eer_projection(:,kk) = circshift(digi_minimum_energy_eer_projection(:,kk),wb_idx(kk));
%     digi_confidence(:,kk) = circshift(digi_confidence(:,kk),wb_idx(kk));
%     if output_std == 1;
%         std_intercept(:,kk) = circshift(std_intercept(:,kk),wb_idx(kk));
%         std_gradient(:,kk) = circshift(std_gradient(:,kk),wb_idx(kk));
%         std_minimum_energy_eer_projection(:,kk) = circshift(std_minimum_energy_eer_projection(:,kk),wb_idx(kk));
%     end
% end

resultno = 1;
% Save outputs into correct structure to be written to SEGY.
results_out{resultno,1} = 'Meta data for output files';
results_out{resultno,2}{1,1} = ilxl_read{vol_index_wb}(1:ntraces,:);
%results_out{resultno,2}{2,1} = uint32(zeros(size(traces{vol_index_wb},2),1));
results_out{resultno,2}{2,1} = uint32(zeros(ntraces,1));

ebcstrtowrite = sprintf('%-3200.3200s',[results_out{resultno,1} '  ' ebdichdr '  ' tmpebc]);
results_out{resultno,1} = ebcstrtowrite;

resultno = resultno + 1;
clear ilxl_read;

results_out{1,3} = 'is_gather'; % 1 is yes, 0 is no

results_out{resultno,1} = strcat('digi_intercept',testdiscpt);
%results_out{2,2} = digi_intercept;
results_out{resultno,2} = zeros(job_meta.n_samples{vol_index_wb},ntraces);
% Unflatten data using the window index
results_out{resultno,2}(win_ind(:,1:ntraces)) = 1000.*ava(1:ns,:);
results_out{resultno,3} = 0;
resultno = resultno + 1;

results_out{resultno,1} = strcat('digi_gradient',testdiscpt);
%results_out{3,2} = digi_gradient;
results_out{resultno,2} = zeros(job_meta.n_samples{vol_index_wb},ntraces);
% Unflatten data using the window index
results_out{resultno,2}(win_ind(:,1:ntraces)) = 1000.*ava(1+ns:end,:);
results_out{resultno,3} = 0;
resultno = resultno + 1;

results_out{resultno,1} = strcat('digi_minimum_energy_eer_projection',testdiscpt); 
%results_out{4,2} = digi_minimum_energy_eer_projection;
digi_minimum_energy_eer_projection = [bsxfun(@times,ava(1:ns,:),cosd(chi))+bsxfun(@times,ava(1+ns:end,:),sind(chi));zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
results_out{resultno,2} = zeros(job_meta.n_samples{vol_index_wb},ntraces);
% Unflatten data using the window index
results_out{resultno,2}(win_ind(:,1:ntraces)) = 1000.*digi_minimum_energy_eer_projection(1:ns,:);
results_out{resultno,3} = 0;
resultno = resultno + 1;

if needconf == 1;
    results_out{resultno,1} = strcat('digi_confidence',testdiscpt);
    %results_out{5,2} = digi_confidence;
    digi_confidence = [digi_confidence;zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
    results_out{resultno,2} = zeros(job_meta.n_samples{vol_index_wb},ntraces);
    % Unflatten data using the window index
    results_out{resultno,2}(win_ind(:,1:ntraces)) = 1000.*digi_confidence(1:ns,:);
    results_out{resultno,3} = 0;
    resultno = resultno + 1;
end

% output standard intercept and gradient result
if output_std == 1;
    results_out{resultno,1} = strcat('std_intercept',testdiscpt);
    results_out{resultno,2} = 1000.*std_intercept;
    results_out{resultno,3} = 0;
    resultno = resultno + 1;
    
    results_out{resultno,1} = strcat('std_gradient',testdiscpt);
    results_out{resultno,2} = 1000.*std_gradient;
    results_out{resultno,3} = 0;
    resultno = resultno + 1;

    results_out{resultno,1} = strcat('std_minimum_energy_eer_projection',testdiscpt);
    results_out{resultno,2} = 1000.*std_minimum_energy_eer_projection;
    results_out{resultno,3} = 0;
    resultno = resultno + 1;
end
%%
% segy write function
if exist(strcat(job_meta.output_dir,'digi_results/'),'dir') == 0
    output_dir = strcat(job_meta.output_dir,'digi_results/');
    mkdir(output_dir);    
else
    output_dir = strcat(job_meta.output_dir,'digi_results/');
end

i_block = str2double(i_block);
node_segy_write(results_out,i_block,job_meta.s_rate/1000,output_dir);

end

%%
function [IGmatrix] = build_operator(totalvol,input_angles,ns_wavelet,wavelet_tmp,ns,hns_wavelet,angle_index,chi,alpha,IGblank,smooth)
    % Build operator for inversion
    for ii = 1:totalvol
        Iwavelet_interp(:,1+(ii-1)*ns_wavelet:ii*ns_wavelet) = wavelet_tmp{ii}';
        %Iwavelet_interp(:,1+(ii-1)*ns_wavelet:ii*ns_wavelet) = wavelet_tmp{ii}'*(cosd(input_angles(ii)).*cosd(input_angles(ii)));
        Gwavelet_interp(:,1+(ii-1)*ns_wavelet:ii*ns_wavelet) = wavelet_tmp{ii}'*(sind(input_angles(ii)).*sind(input_angles(ii)));
        %Iwavelet_interp(:,1+(ii-1)*ns_wavelet:ii*ns_wavelet) = interp1(wavelet_z_grid,wavelet_tmp{ii}',start_interp:1:end_interp,'linear','extrap');
        %Gwavelet_interp(:,1+(ii-1)*ns_wavelet:ii*ns_wavelet) = interp1(wavelet_z_grid,...
        %    wavelet_tmp{ii}'*(sind(input_angles(ii)).*sind(input_angles(ii))),start_interp:1:end_interp,'linear','extrap');
    end

    IGdiagnals = sort(reshape([(-hns_wavelet:hns_wavelet)',bsxfun(@plus,(-hns_wavelet:hns_wavelet)',(-ns:-ns:-ns*(totalvol-1)))],1,[]),'descend');

    Imatrix = spdiags(Iwavelet_interp,IGdiagnals,ns*totalvol,ns);
    Gmatrix = spdiags(Gwavelet_interp,IGdiagnals,ns*totalvol,ns);

    EERmatrix = alpha*[bsxfun(@times,Imatrix(1+ns*(angle_index-1):+ns*angle_index,:),cosd(chi)),bsxfun(@times,Imatrix(1+ns*(angle_index-1):+ns*angle_index,:),sind(chi))];

    IGmatrix = [[Imatrix,Gmatrix;EERmatrix].*IGblank;smooth];
end
% 
% function wavelet_interp = interpolate_wavelets(wavelet_z_grid,wavelet,start_interp,end_interp)
% %     start_interp = min(wavelet_z_grid)-mode(diff(wavelet_z_grid'));
% %     end_interp = max(wavelet_z_grid)+mode(diff(wavelet_z_grid'));
%     wavelet_interp = interp1(wavelet_z_grid,wavelet',start_interp:1:end_interp,'spline','extrap');
% end
%%
function wavelet_norm = wavelet_rms_norm(totalvol,wavelet_interp,vol_index)
    % Normalise the wavelets to have constant energy w.r.t. angle. The energy
    % is set to that of the nearest angle wavelets. Wavelet energy still varies
    % w.r.t. time.
    norm_to = sum(abs(wavelet_interp{vol_index}));
    for ii=1:totalvol
        curwav = sum(abs(wavelet_interp{ii}));
        ratio = norm_to./curwav;
        wavelet_norm{ii} = bsxfun(@times,wavelet_interp{ii},ratio);
    end
%     A = cell2mat(wavelet_interp);
%     %B = sqrt(sum(A.^2));
%     B = (sum(abs(A)));
%     %B = sqrt(max(A.^2));
%     C = reshape(B',length(wavelet_z_grid),[]);
%     D = C(:,vol_index);
%     for ii=1:totalvol
%         E = A(:,1+(ii-1)*length(wavelet_z_grid):ii*length(wavelet_z_grid));
%         F = bsxfun(@rdivide,bsxfun(@times,E,D'),sqrt(sum(E.^2)));
%         F(isnan(F)) = 0;
%         wavelet_norm{ii} = F;
%     end
end
