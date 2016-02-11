function [] = residual_velocity_analysis_aruba( gather_meta_path,i_block,outdir,vel_meta_path,varargin)
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

% Automatic residual velocity analysis
%   
% Picks residual time errors on NMO corrected gathers to give updated
% velocity field
%
% gather_meta_path:     pathname to .mat file for the gathers
% i_block:              block number to analyse
% outdir:               directory for output files
% vel_meta:             path to velocity meta file
%                       note vels must be on same il,xl grid as gathers
%
% Optional parameters:
%
% nmo=                  should be set to 1 if input gathers do not have nmo
%                       applied (default 0)
% gap=                  min gap between picks (number of samples) 
%                       (default 4)
% itm=                  inner trace mute (degrees) (default 0)
% otm=                  outer trace mute (degrees) (default 30)
% min_fold=             minimum number of near offset traces to keep live
%                       (default 5)
% threshold=            minimum cross-correlation value for picks (0 - 1)
%                       (default 0.9)
% start=                horizon at which to start analysis 
%                       (default start at zero)
% stop=                 horizon at which to stop analysis
%                       (should be added with add_horizon_to_job_meta)
%                       (default end at end of trace)
%
% velocity field should be sampled on same x,y grid as gathers but can be
% decimated in t
%
% all optional parameters can be left to default to get a reasonable result
% but adding start and stop horizons can greatly reduce runtime

% to reduce printout in compilied version turned all warning off
warning off all;

% Defaults for arguments =================================================
apply_nmo_switch = 0;
gap = 4;
itm = 0;
otm = 30;
min_fold = 5;
threshold = 0.9;
use_start_horz = 0;
start_horz = 'zero';
use_stop_horz = 0;
stop_horz = 'eot';
i_block = str2num(i_block);
%=========================================================================
%
for kv = 1:length(varargin)
    varknown = false;
    if strfind(varargin{kv},'nmo=')
        vartmp = deblank(regexp(varargin{kv},'=','split'));
        apply_nmo_switch = str2double(vartmp(2));
        varknown = true;
    end
    if strfind(varargin{kv},'gap=')
        vartmp = deblank(regexp(varargin{kv},'=','split'));
        gap = str2double(vartmp(2));
        varknown = true;
    end
    if strfind(varargin{kv},'itm=')
        vartmp = deblank(regexp(varargin{kv},'=','split'));
        itm = str2double(vartmp(2));
        varknown = true;
    end
    if strfind(varargin{kv},'otm=')
        vartmp = deblank(regexp(varargin{kv},'=','split'));
        otm = str2double(vartmp(2));
        varknown = true;
    end
    if strfind(varargin{kv},'min_fold=')
        vartmp = deblank(regexp(varargin{kv},'=','split'));
        min_fold = str2double(vartmp(2));
        varknown = true;
    end
    if strfind(varargin{kv},'threshold=')
        vartmp = deblank(regexp(varargin{kv},'=','split'));
        threshold = str2double(vartmp(2));
        varknown = true;
    end
    if strfind(varargin{kv},'start=')
        vartmp = deblank(regexp(varargin{kv},'=','split'));
        use_start_horz = 1;
        start_horz = vartmp{2};
        varknown = true;
    end
    if strfind(varargin{kv},'stop=')
        vartmp = deblank(regexp(varargin{kv},'=','split'));
        use_stop_horz = 1;
        stop_horz = vartmp{2};
        varknown = true;
    end
    
    % last test to see if there was a command line variable that was not
    % expected in which case throw warning and crash
    if varknown == false;
        error('unknown input variable; the variable %s is not recognised and the function has quit\n', varargin{kv});
    end
end


% ============================================================================================================
%
% load the data
%

gather_meta = load(gather_meta_path);
vel_meta = load(vel_meta_path);

min_il = gather_meta.block_keys(i_block,1);
max_il = gather_meta.block_keys(i_block,2);
min_xl = gather_meta.block_keys(i_block,3);
max_xl = gather_meta.block_keys(i_block,4);

il_inc = double(gather_meta.pkey_inc);
xl_inc = double(gather_meta.skey_inc);

num_ils = 1 + (max_il-min_il)/il_inc;
num_xls = 1 + (max_xl-min_xl)/xl_inc;

% do some basic checks

error_flag = 0;
error_str = '';
if vel_meta.pkey_min > min_il
    error_str = ['Velocity volume first inline (',num2str(vel_meta.pkey_min),') greater than first inline for block ',num2str(i_block),' (',num2str(min_il),')'];
    error_flag = 1;
end
if vel_meta.pkey_max < max_il
    error_str = [error_str,'Velocity volume last inline (',num2str(vel_meta.pkey_max),') less than last inline for block ',num2str(i_block),' (',num2str(max_il),')'];
    error_flag = 1;
end
if vel_meta.skey_min > min_xl
    error_str = [error_str,'Velocity volume first xline (',num2str(vel_meta.skey_min),') greater than first xline for block ',num2str(i_block),' (',num2str(min_xl),')'];
    error_flag = 1;
end
if vel_meta.skey_max < max_xl
    error_str = [error_str,'Velocity volume last xline (',num2str(vel_meta.skey_max),') less than last xline for block ',num2str(i_block),' (',num2str(max_xl),')'];
    error_flag = 1;
end
if vel_meta.pkey_inc ~= il_inc
    error_str = [error_str,'Velocity volume inline increment (',num2str(vel_meta.pkey_inc),') not equal to gather inline increment ' (',num2str(il_inc),')'];
    error_flag = 1;
end
if vel_meta.skey_inc ~= xl_inc
    error_str = [error_str,'Velocity volume xline increment (',num2str(vel_meta.skey_inc),') not equal to gather xline increment ' (',num2str(xl_inc),')'];
    error_flag = 1;
end

if error_flag
    error(error_str);
end



% ============================================================================================================
%
% read the velocity volume
%
% find the blocks to read
counter=0;

for bb=1:length(vel_meta.liveblocks);
    vel_blk=vel_meta.liveblocks(bb);
    vel_min_il=vel_meta.block_keys(vel_blk,1);
    vel_max_il=vel_meta.block_keys(vel_blk,2);
    vel_min_xl=vel_meta.block_keys(vel_blk,3);
    vel_max_xl=vel_meta.block_keys(vel_blk,4);
    vel_il_corners = [vel_min_il,vel_min_il,vel_max_il,vel_max_il];
    vel_xl_corners = [vel_min_xl,vel_max_xl,vel_max_xl,vel_min_xl];
    in_block = inpolygon(vel_il_corners,vel_xl_corners,[min_il,min_il,max_il,max_il],[min_xl,max_xl,max_xl,min_xl]);
    if sum(in_block)>0
        counter=counter+1;
       live_velblocks(counter)=vel_blk;
    end
end

veltrace = zeros(vel_meta.n_samples{1},num_ils*num_xls);

for bb=1:counter
    [velstr, veltrace_blk, veltrace_ilxl_blk, ~] = node_segy_read(vel_meta_path,'1',num2str(live_velblocks(bb)));
    
    keep_idx = and(and(veltrace_ilxl_blk(:,1) >= min_il,veltrace_ilxl_blk(:,1) <= max_il),and(veltrace_ilxl_blk(:,2) >= min_xl,veltrace_ilxl_blk(:,2) <= max_xl));
    
    veltrace_keep = veltrace_blk(:,keep_idx);
    veltrace_ilxl_keep = veltrace_ilxl_blk(keep_idx,:);
    
    il_idx = int32((veltrace_ilxl_keep(:,1)-min_il)./il_inc);
    xl_idx = int32((veltrace_ilxl_keep(:,2)-min_xl)./xl_inc);
    ilxl_idx = il_idx.*num_xls + xl_idx + 1;
    
    veltrace(:,ilxl_idx') = veltrace_keep;
end

clear veltrace_blk veltrace_keep veltrace_ilxl_blk veltrace_ilxl_keep il_idx xl_idx keep_idx

veltimes = 0.001*((velstr.s_rate/1000):(velstr.s_rate/1000):(velstr.s_rate/1000)*velstr.n_samples)';

% read the gathers and velocities in from segy

% trace data is ordered by cdp, offset, sample

% incoming trace data is one stream of traces for efficiency
% so reshape it into a 3 dimensional array for processing

[seismic, traces, traces_ilxl, offsets] = node_segy_read(gather_meta_path,'1',num2str(i_block));

seismic.fold = max(seismic.trace_ilxl_bytes(:,7));
traces = reshape(traces,size(traces,1),seismic.fold,[]);


seismic.n_gathers = size(traces,3);
% check to make sure it read something if not exit

if isempty(traces) == 1 &&  isempty(traces_ilxl) == 1 && isempty(offsets) == 1
    return
end

gather_ilxl = unique(traces_ilxl,'rows');

offsets = reshape(offsets,seismic.fold,[]);

% ============================================================================================================
% load start and stop horizons

% initialise start and stop times
horz_times = zeros(seismic.n_gathers,2);
horz_times(:,2) = seismic.n_samples*seismic.s_rate/1000;

if use_start_horz
    horz_in = dlmread(gather_meta.(start_horz));
    keep_idx = and(and(horz_in(:,1) >= min_il,horz_in(:,1) <= max_il),and(horz_in(:,2) >= min_xl,horz_in(:,2) <= max_xl));
    horz_keep = horz_in(keep_idx,:);
    il_idx = int32((horz_keep(:,1)-min_il)./il_inc);
    xl_idx = int32((horz_keep(:,2)-min_xl)./xl_inc);
    ilxl_idx = il_idx.*num_xls + xl_idx + 1;
    horz_times(ilxl_idx,1) = horz_keep(:,3);
end

clear horz_in keep_idx horz_keep il_idx xl_idx ilxl_idx;

if use_stop_horz
    horz_in = dlmread(gather_meta.(stop_horz));
    keep_idx = and(and(horz_in(:,1) >= min_il,horz_in(:,1) <= max_il),and(horz_in(:,2) >= min_xl,horz_in(:,2) <= max_xl));
    horz_keep = horz_in(keep_idx,:);
    il_idx = int32((horz_keep(:,1)-min_il)./il_inc);
    xl_idx = int32((horz_keep(:,2)-min_xl)./xl_inc);
    ilxl_idx = il_idx.*num_xls + xl_idx + 1;
    horz_times(ilxl_idx,2) = horz_keep(:,3);
end


% % ============================================================================================================
% 
% % set up meta data for output
% 
% % add the history of jobs run and this one to the current ebcdic
% 
% ebdichdr = [strcat('residual velocity analysis ',date)];
% if isfield(gather_meta,'comm_history')
%     prev_ebcdic = gather_meta.comm_history;
%     new_ebcdic = prev_ebcdic{size(prev_ebcdic,1),2};
% else
%     prev_ebcdic{1,2} = '';
%     new_ebcdic = '';
% end
% 
% for ebcii = (size(prev_ebcdic,1)-1):-1:1
%     tmpebcc = regexp(prev_ebcdic{ebcii,2},'/','split');
%     new_ebcdic = [new_ebcdic tmpebcc{1}  tmpebcc{end}];
% end
% new_ebcdic = sprintf('%-3200.3200s',new_ebcdic);
% clear tmpebcc prev_ebcdic;
% 
% % for the stack output
% stack_traces{1,1} = 'Meta data for output files';
% stack_traces{1,2}{1,1} = traces_ilxl(1:seismic.fold:end,:);
% stack_traces{1,2}{2,1} = uint32(zeros(seismic.n_gathers,1));
% %ebcstrtowrite = sprintf('%-3200.3200s',[traces{resultno,1} '  ' ebdichdr '  ' tmpebc]);
% stack_traces{1,1} = ebcstrtowrite;
% stack_traces{1,3} = 'is_gather'; % 1 is yes, 0 is no
% 
% output_dir = [gather_meta.output_dir,'trimout/'];
% % check to see if the directory exists
% if exist(output_dir,'dir') == 0
%     mkdir(output_dir);
% end
% 
% 
% 
% % stack result
% filename2 = 'stack';
% stack_traces{2,1} = strcat(filename2,'_input');
% stack_traces{2,2} = zeros(seismic.n_samples,seismic.n_gathers,'single');
% stack_traces{3,1} = strcat(filename2,'_output');
% stack_traces{3,2} = zeros(seismic.n_samples,seismic.n_gathers,'single');
% %stack_traces{4,1} = strcat(filename2,'_posttrim_no_resid_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
% stack_traces{4,1} = strcat(filename2,'_velocity');
% stack_traces{4,2} = zeros(seismic.n_samples,seismic.n_gathers,'single');
% stack_traces{2,3} = 0;
% stack_traces{3,3} = 0;
% stack_traces{4,3} = 0;

% ============================================================================================================
% apply nmo

if apply_nmo_switch == 1
   for gather = 1:seismic.n_gathers
       traces(:,:,gather) = nmo(traces(:,:,gather),seismic.s_rate/1000000,offsets(:,gather),veltimes(:,gather),veltrace(:,gather),1000);
   end
end

 

% ============================================================================================================
% calculate offsets from angles and velocities and make mute masks

angles = [itm otm];

nsamples = size(veltrace,1);

pick_mute=zeros(nsamples,seismic.fold,seismic.n_gathers);

ntraces=seismic.n_gathers;

offset_traces = zeros(nsamples,size(angles),ntraces);

angles_rad = 2*pi*angles/360;

% =====================================================================

% angle calculation


% tan theta = 0.5 * offset / depth = 0.5 * offset / (0.5 * twtt * vrms)
% => offset = tan theta * twtt * vrms

for ang_idx=1:max(size(angles))
    
    offset_traces(:,ang_idx,:) = bsxfun(@times,veltimes,veltrace * tan(angles_rad(ang_idx)));
    
end

% compare trace offsets with those calculated for the mute angles and
% create a mute mask for each gather

itm_mute = squeeze(offset_traces(:,1,:));
otm_mute = squeeze(offset_traces(:,2,:));

% make matrix to keep n traces live when applying otm

otm_limit = [ones(min_fold,1);zeros(seismic.fold-min_fold,1)];
   

for i_samp = 1:size(offset_traces,1)    
    mute_mask(i_samp,:,:) = bsxfun(@gt,offsets,itm_mute(i_samp,:)).*bsxfun(@or,bsxfun(@lt,offsets,otm_mute(i_samp,:)),otm_limit);
end

% find range of live traces

ssr = seismic.s_rate/1000;
sn = seismic.n_samples;

vs=velstr.s_rate/1000;
vn=nsamples;

% interp mute mask onto seismic sampling
% pad mute with extra row at start and interpolate:

mute_mask = interp1([ssr vs:vs:vs*vn vs*vn+vs],[mute_mask(1,:,:);mute_mask;mute_mask(end,:,:)],ssr:ssr:ssr*sn);

% pad and smooth mutes then remove padding

taperlen = 1+40/ssr;
smth = linspace(0,1,0.5*taperlen);
smth = [smth smth((end-1):-1:2)];
smth = smth'/sum(smth);

mute_mask = [mute_mask;repmat(mute_mask(end,:,:),size(smth,1)*2,1)];
mute_mask = convn(mute_mask,smth,'same');
mute_mask = mute_mask(1:seismic.n_samples,:,:);

% apply the mute
traces = traces.*mute_mask;

clear ssr; clear sn; clear vs; clear vn;

% ============================================================================================================
% pick trim shifts and get maximum x-correlation peaks

% we want output to be inline,xline,time,velocity
gap = 12;
output_matrix = zeros((seismic.n_samples/gap)*seismic.n_gathers,4);
outcount = 0;
% run just for a subset for testing
for gather = 1:seismic.n_gathers
    
    gather_in = traces(:,:,gather);
    start_idx = round(1000*horz_times(gather,1)./seismic.s_rate);
    stop_idx = round(1000*horz_times(gather,2)./seismic.s_rate);
    vert_sum = sum(gather_in(start_idx:stop_idx,:),1);
    start_gath = find(vert_sum,1,'first');
    end_gath = find(vert_sum,1,'last');
    
    gather_in = gather_in(start_idx:stop_idx,start_gath:end_gath);
        
    [trim_shifts,~,~,trim_coeffs] = trim_shifts_calculate('calc',gather_meta,gather_in,'1','4');
    
    % pick peak values of xcor coefficients to get times for picks
    %
    % this loop picks the maximum value, then sets surrounding values
    % to zero and continues to loop until there are no more values over
    % the minimum threshold. Might be faster way to do it.
    
    coeff_peaks_out = zeros(size(trim_coeffs));
    
    while max(trim_coeffs)>threshold % threshold is a user parameter
        [~,peak_idx] = max(trim_coeffs);
        coeff_peaks_out(peak_idx) = trim_coeffs(peak_idx);
        start_mask=max(peak_idx-gap,1);
        end_mask=min(peak_idx+gap,seismic.n_samples);
        trim_coeffs(start_mask:end_mask) = 0;
    end
    
    pick_idx = find(coeff_peaks_out);
    
    ptimes = (start_idx + pick_idx) * seismic.s_rate/1000;
    
    %==============================================================
    
    numpicks = size(ptimes,1);
    outcount = outcount + numpicks;
    vels = 0.001*interp1(veltimes*1000,veltrace(:,gather),ptimes); % times are in ms, so make vels m/ms
    
    t0times = bsxfun(@minus,trim_shifts(pick_idx,:),ptimes);
    ttimes = sqrt(t0times.^2 + bsxfun(@rdivide,double(offsets(start_gath:end_gath,gather)'.^2),vels.^2));
    
    % set times in mute zone to be NaN
    
    ttimes(mute_mask(pick_idx+start_idx,start_gath:end_gath,gather)<1) = NaN;
    vel_out = zeros(numpicks,1);
    for ii = 1:numpicks
        tsq = (ttimes(ii,~isnan(ttimes(ii,:)))./1000).^2;
        offsq = offsets(~isnan(ttimes(ii,:)),1000).^2;
        p = polyfit(double(offsq),tsq',1);
        vel_out(ii) = 1/sqrt(p(1));
        %             vel_out(gather,ii) = 1/sqrt(p(1));
        
    end
    
    output_matrix(outcount-numpicks+1:outcount,1) = gather_ilxl(gather,1);
    output_matrix(outcount-numpicks+1:outcount,2) = gather_ilxl(gather,2);
    output_matrix(outcount-numpicks+1:outcount,3) = ptimes;
    output_matrix(outcount-numpicks+1:outcount,4) = vel_out;
       
end

% vint = bsxfun(@velfun,veltrace,veltimes);

out_filename = strcat(outdir,'rms_picks_gap',num2str(gap),'_mute',num2str(itm),'_',num2str(otm),'_thr',num2str(threshold*100),'_',start_horz,'_',stop_horz,'_block',num2str(i_block));

dlmwrite(strcat(out_filename,'.txt'),output_matrix(1:outcount,:));

save(strcat(out_filename,'.mat'),'output_matrix');

% write some output


% stack_traces{2,2} = stack;
% stack_traces{4,2} = velout;
% i_block = str2double(i_block);
% node_segy_write(stack_traces,i_block,seismic.s_rate/1000,outdir);
% 
% save(strcat(outdir,'rms_picks_block',int2str(i_block),'.mat'),'ntraces','tvpairs','tvpairs2',...
%     'stack_traces','stack','stack_mask','resdiff','resdiffnorm','resdiff_signed','resdiff_signed_norm');


end

function vint = velfun(vrms,ttimes)
tts = [0;ttimes(1:end-1)];
vrms_s = [vrms(1);vrms(1:end-1)];

vint = sqrt(((vrms.^2.*ttimes) - (vrms_s.^2.*tts)) ./ (ttimes-tts));
end

