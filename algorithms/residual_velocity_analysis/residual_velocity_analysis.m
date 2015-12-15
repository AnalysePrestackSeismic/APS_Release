function [] = residual_velocity_analysis( gather_meta_path,i_block,outdir,apply_nmo_switch)

% Automatic residual velocity analysis
%   
% Inverts residual time errors on NMO corrected gathers to give updated
% velocity field
%
% gather_meta_path:     pathname to .mat file for the gathers
% i_block:              block number to analyse
% outdir:               directory for output files
% apply_nmo_switch:     should be set to 1 if need to apply nmo, otherwise 0           
%
% velocity field should be sampled on same x,y grid as gathers but can be
% decimated in t
%

% to reduce printout in compilied version turned all warning off
warning off all;
apply_nmo_switch = str2double(apply_nmo_switch);

% ============================================================================================================
% hardwired parameters

stack_itm = 5;
stack_otm = 80;
pick_otm = 70;

% num_iter=5;
% smth_size=9;


% ============================================================================================================
%
% load the data
%

gather_meta = load(gather_meta_path);

% make a single velocity trace for testing with (average t/s dip profile)


% veltrace = interp1([0 132 264 396 529 663 797 1068 1607 2676 3341 4005 4667 5330 5992 6654],...
%     [1525 1520 1517 1515 1512 1509 1505 1499 1493 1495 1497 1498 1500 1501 1502 1503],100:100:5000);

% make a single velocity trace for testing with (constant 1500 m/s)

veltrace = interp1([0 6000],[1500 1500],100:100:5000);

veltimes = ([100:100:5000]')./1000;

veltrace = veltrace';


% velocity_meta = load(velocity_meta);

% read the gathers and velocities in from segy

% trace data is ordered by cdp, offset, sample

% incoming trace data is one stream of traces for efficiency
% so reshape it into a 3 dimensional array for processing

[seismic, traces{3,2}, traces{1,2}{1,1}, offsets] = node_segy_read(gather_meta_path,'1',i_block);


seismic.fold = max(seismic.trace_ilxl_bytes(:,7));
traces{3,2} = reshape(traces{3,2},size(traces{3,2},1),seismic.fold,[]);

% [velocity, vel_traces{3,2}, vel_traces{1,2}{1,1},~] = node_segy_read(velocity_meta_path,'0',i_block)

seismic.n_gathers = size(traces{3,2},3);
% check to make sure it read something if not exit

veltrace = repmat(veltrace,1,seismic.n_gathers);
veltimes = repmat(veltimes,1,seismic.n_gathers);


if isempty(traces{3,2}) == 1 &&  isempty(traces{1,2}{1,1}) == 1 && isempty(offsets) == 1
    return
end
%if isempty(vel_traces{3,2}) == 1 &&  isempty(vel_traces{1,2}{1,1}) == 1 
%    return
%end

% ============================================================================================================

% set up meta data for output

% add the history of jobs run and this one to the current ebcdic

ebdichdr = [strcat('residual velocity analysis ',date)];
if isfield(gather_meta,'comm_history')
    prev_ebcdic = gather_meta.comm_history;
    new_ebcdic = prev_ebcdic{size(prev_ebcdic,1),2};
else
    prev_ebcdic{1,2} = '';
    new_ebcdic = '';
end

for ebcii = (size(prev_ebcdic,1)-1):-1:1
    tmpebcc = regexp(prev_ebcdic{ebcii,2},'/','split');
    new_ebcdic = [new_ebcdic tmpebcc{1}  tmpebcc{end}];
end
new_ebcdic = sprintf('%-3200.3200s',new_ebcdic);
clear tmpebcc prev_ebcdic;


% for the pre-stack dataset
traces{1,1} = 'Meta data for output files';
%traces{resultno,2}{1,1} = ilxl_read;
traces{1,2}{2,1} = offsets;
ebcstrtowrite = sprintf('%-3200.3200s',[traces{1,1} '  ' ebdichdr '  ' new_ebcdic]);
traces{1,1} = ebcstrtowrite;
traces{1,3} = 'is_gather'; % 1 is yes, 0 is no


% for the stack output
stack_traces{1,1} = 'Meta data for output files';
stack_traces{1,2}{1,1} = traces{1,2}{1,1}(1:seismic.fold:end,:);
stack_traces{1,2}{2,1} = uint32(zeros(seismic.n_gathers,1));
%ebcstrtowrite = sprintf('%-3200.3200s',[traces{resultno,1} '  ' ebdichdr '  ' tmpebc]);
stack_traces{1,1} = ebcstrtowrite;
stack_traces{1,3} = 'is_gather'; % 1 is yes, 0 is no

output_dir = [gather_meta.output_dir,'trimout/'];
% check to see if the directory exists
if exist(output_dir,'dir') == 0
    mkdir(output_dir);
end


% prestack output
filename = 'gaths';
%traces{2,1} = strcat(filename,'_trim_shifts_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
traces{2,1} = strcat(filename,'_trim_shifts');
traces{3,1} = strcat(filename,'_nmo_data');
traces{2,3} = 1;
traces{3,3} = 1;
traces{2,2} = zeros(size(traces{3,2}),'single');


% stack result
filename2 = 'stack';
stack_traces{2,1} = strcat(filename2,'_input');
stack_traces{2,2} = zeros(seismic.n_samples,seismic.n_gathers,'single');
stack_traces{3,1} = strcat(filename2,'_output');
stack_traces{3,2} = zeros(seismic.n_samples,seismic.n_gathers,'single');
%stack_traces{4,1} = strcat(filename2,'_posttrim_no_resid_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
stack_traces{4,1} = strcat(filename2,'_velocity');
stack_traces{4,2} = zeros(seismic.n_samples,seismic.n_gathers,'single');
stack_traces{2,3} = 0;
stack_traces{3,3} = 0;
stack_traces{4,3} = 0;

% ============================================================================================================
% apply nmo

offsets = reshape(offsets,seismic.fold,[]);

if apply_nmo_switch == 1
   for gather = 1:seismic.n_gathers
       traces{3,2}(:,:,gather) = nmo(traces{3,2}(:,:,gather),seismic.s_rate/1000000,offsets(:,gather),veltimes(:,gather),veltrace(:,gather),1000);
   end
end

 

% ============================================================================================================
% calculate offsets from angles and velocities and make mute masks

angles = [stack_itm stack_otm pick_otm];

nsamples = size(veltrace,1);

stack_mute=zeros(nsamples,seismic.fold,seismic.n_gathers);
pick_mute=zeros(size(stack_mute));

% stack_mute=zeros(size(vel_traces{3,2}));
% pick_mute=zeros(size(vel_traces{3,2}));
% stack_fold=zeros(velocity.nsamples,velocity.n_gathers);
% pick_fold=zeros(velocity.nsamples,velocity.n_gathers);

ntraces=seismic.n_gathers;

offset_traces = zeros(nsamples,size(angles),ntraces);

angles_rad = 2*pi*angles/360;

% ======================================================================
%
% dix convert rms to interval
%
% replace with constrained inversion
%
vint=zeros(nsamples,ntraces);
tt = ((1:nsamples).*100)';

%tts = [0;tt(1:end-1)];
%veltrace_s = [veltrace(1);veltrace(1:end-1)];

%vint = sqrt(((veltrace.^2.*tt) - (veltrace_s.^2.*tts)) ./ (tt-tts));

% this is how to do it if you have more than one vel trace but only one
% column of times

vint = bsxfun(@velfun,veltrace,tt);



% vint = sqrt(bsxfun(@rdivide,(bsxfun(@times,(veltrace.^2),tt) - bsxfun(@times,(veltrace_s.^2),tts)),(tt-tts)));



% =====================================================================

% angle calculation

% sine theta = (offset * vint) / (t * vrms^2)

% => offset = (t * vrms^2 * sine theta) / vint


for ang_idx=1:max(size(angles))
    
    offset_traces(:,ang_idx,:) = bsxfun(@times,tt,((veltrace/1000).^2 * sin(angles_rad(ang_idx)))) ./ (vint/1000);
    
%     for trc_idx=1:ntraces
%     
%         offset_traces(:,ang_idx,trc_idx) = (tt(:) .* vel_traces(:,trc_idx) .* sin(angles_rad(ang_idx))) ./ vint(:,trc_idx);
%     
%     end
end


% compare trace offsets with those calculated for the mute angles and
% create a mute mask for each gather

mute1 = squeeze(offset_traces(:,1,:));
mute2 = squeeze(offset_traces(:,2,:));
mute3 = squeeze(offset_traces(:,3,:));

% make matrix to keep n traces live when applying otm

otm_limit = [ones(7,1);zeros(seismic.fold-7,1)];
   

for i_samp = 1:1:size(offset_traces,1)    
    stack_mute(i_samp,:,:) = bsxfun(@gt,offsets,mute1(i_samp,:)).*bsxfun(@or,bsxfun(@lt,offsets,mute2(i_samp,:)),otm_limit);
    pick_mute(i_samp,:,:) = bsxfun(@lt,offsets,mute3(i_samp,:));
end


% for gather_idx=1:seismic.ngathers
%     
%     % offset_lookup = offset_vs_angle(velocity_meta,angles,vel_traces{3,2});
%     
%     for off_idx=1:seismic.fold
%         stack_mute(:,off_idx,gather_idx) = (offsets(gather_idx) >= offset_lookup(:,1,gather_idx)) ...
%             .* (offsets(gather_idx) <= offset_lookup(:,2,gather_idx));
%         pick_mute(:,off_idx,gather_idx) = (offsets(gather_idx) <= offset_lookup(:,3,gather_idx));
%     end
%     stack_fold(:,gather_idx) = sum(stack_mute,2);
%     pick_fold(:,gather_idx) = sum(pick_mute,2);
% end

% interpolate the mutes onto the seismic sampling rate

ssr = seismic.s_rate/1000;
sn = seismic.n_samples;
%vs = velocity.srate;
%vn = velocity.nsamples;

vs=100;
vn=nsamples;

stack_mute = interp1([ssr vs:vs:vs*vn vs*vn+vs],[stack_mute(1,:,:);stack_mute;stack_mute(end,:,:)],ssr:ssr:ssr*sn);
pick_mute = interp1([ssr vs:vs:vs*vn vs*vn+vs],[pick_mute(1,:,:);pick_mute;pick_mute(end,:,:)],ssr:ssr:ssr*sn);


% interp mute mask onto seismic sampling


% pad mute with extra row at start and interpolate:


% stack_traces{2,2} = (sum(traces{3,2},2);

% pad and smooth mutes then remove padding

taperlen = 1+40/ssr;
smth = linspace(0,1,0.5*taperlen);
smth = [smth smth((end-1):-1:2)];
smth = smth'/sum(smth);

stack_mute = [stack_mute;repmat(stack_mute(end,:,:),size(smth,1)*2,1)];
stack_mute = convn(stack_mute,smth,'same');
stack_mute = stack_mute(1:seismic.n_samples,:,:);

pick_mute = [pick_mute;repmat(pick_mute(end,:,:),size(smth,1)*2,1)];
pick_mute = convn(pick_mute,smth,'same');
pick_mute = pick_mute(1:seismic.n_samples,:,:);


clear ssr; clear sn; clear vs; clear vn;

% ============================================================================================================
% pick trim shifts and get maximum x-correlation peaks

unmuted_traces = traces{3,2};

traces{3,2} = traces{3,2}.*stack_mute;

% run just for a subset for testing

[trim_shifts,~,~,trim_coeffs] = trim_shifts_calculate('calc',gather_meta,traces{3,2},'1','4');

% pick xcorrelation peaks - these are times to pick velocities at

% fairly sure this is v inefficient way of doing it
% also try smoothing the coefficents to get more consistent time picks

% mask out areas for picking based on the stack
stack = squeeze(sum(traces{3,2},2));
stackfold = squeeze(sum(stack_mute,2));
stack = stack ./ stackfold;
% gaussian
t = linspace(-1,1,51)';
filttraces = zeros(length(t),1);
a = 1;
filttraces = sqrt(pi)/a*exp(-(pi*t/a).^2);
filttraces = filttraces/sum(filttraces);
filttraces = filttraces';

stack_mask=zeros(size(trim_coeffs));
stack_mask(1:2500,:)=bwareaopen(convn(stk_event_zones(stack,10,7),filttraces,'same')>0.4,400);


trim_coeffs = trim_coeffs.*stack_mask;

% zero the first 100 samples

trim_coeffs(1:100,:) = 0;

gap = 20;

coeff_peaks_out = zeros(size(trim_coeffs));

smth = [1];
smth = smth / sum(smth);
pad = ceil(size(smth,2)/2);

coeffs_smth = convn([repmat(trim_coeffs(:,1),1,pad) trim_coeffs repmat(trim_coeffs(:,size(trim_coeffs,2)),1,pad)],smth,'same');
coeffs_smth = coeffs_smth(:,(pad+1):(size(coeffs_smth,2)-pad));


for trc = 1:size(coeffs_smth,2);
    
    while max(coeffs_smth(:,trc))>0.8
        
        [~,peak_idx] = max(coeffs_smth(:,trc));
        coeff_sm_peaks_out(peak_idx,trc) = coeffs_smth(peak_idx,trc);
        start_mask=max(peak_idx-gap,1);
        end_mask=min(peak_idx+gap,seismic.n_samples);
        
        coeffs_smth(start_mask:end_mask,trc) = 0;
        
    end
    
end


% ============================================================================================================

% calc travel times for trim picks (assume straight ray)

% tt = sqrt(t0^2 + (offset^2/v^2))


velout = zeros(1500,seismic.n_gathers,3);


for gather = 1:seismic.n_gathers
%for gather = 1:100
    
  pick_idx = find(coeff_sm_peaks_out(:,gather));
  pick_times = pick_idx * seismic.s_rate/1000;
  
  vels = 0.001*interp1([0 100:100:5000],[veltrace(1,gather);veltrace(:,gather)]',pick_times);
  
  t0times = bsxfun(@minus,trim_shifts(pick_idx,:,gather),pick_times);
  %t0times2 = bsxfun(@plus,-trim_shifts(pick_idx,:,gather),pick_times);
  %t0times2 = bsxfun(@minus,trim_shifts(pick_idx,:,gather),pick_times);
  ttimes = sqrt(t0times.^2 + bsxfun(@rdivide,double(offsets(:,gather)'.^2),vels.^2));

  % set times in mute zone to be NaN
  
  ttimes(stack_mute(pick_idx,:,1)<1) = NaN;
  
  % scale the travel times by the number of picks on each event -
  % this has the effect of weighting each event equally in the solution
  % also scale by 1/t so that the error is proportional to the velocity
  % error
  %
  % also supply initial guess at solution so lsqr has less work to do and
  % is less likely to stagnate before getting to sensible answer
  
  scale_fn = sum(isfinite(ttimes),2);
  scale_fn = scale_fn./pick_times;
  scale_fn = max(scale_fn)./scale_fn;
  ttimes = bsxfun(@times,ttimes,scale_fn);
  ini_tv = [pick_times ones(size(pick_times),1).*1.5];
  ini_tv(:,1) = ini_tv(:,1) .* scale_fn;
  ini_tv(:,2) = ini_tv(:,2) ./ scale_fn;
  
%   smth = [1:ceil(smth_size/2) floor(smth_size/2):-1:1];
%   smth = smth'./sum(smth);
%   pad = floor(size(smth,1)/2);
  
%   for iter = 1:num_iter
%       
      [v_est{gather} t_est{gather}] = invert_T_V_from_picks(ttimes,offsets(:,gather),ini_tv);
      
      v_est{gather} = v_est{gather}.*scale_fn;
      t_est{gather} = t_est{gather}./scale_fn;
      
%       slness = v_est{gather}.^-1;
%       slness = [repmat(slness(1),pad,1);slness;repmat(slness(end),pad,1)];
%       slness = conv(slness,smth,'valid');
%       ini_tv(:,2) = (slness.^-1) ./ scale_fn;

      t_est{gather}(1)=t_est{gather}(1)+1;
      velout(:,gather,1) = interp1([0;t_est{gather}],[v_est{gather}(1);v_est{gather}],2:2:3000);
      
      outtimes = sqrt(bsxfun(@plus,t_est{gather}.^2,bsxfun(@rdivide,double(offsets(:,gather)'.^2),v_est{gather}.^2)));
      intimes = bsxfun(@rdivide,ttimes,scale_fn);
      timediff=outtimes-intimes;
      timediff(isnan(ttimes))=NaN;
      resdiff{gather} = nanmean(abs(timediff),2);
      resdiff_signed{gather} = nanmean(timediff,2);
      resdiffnorm{gather} = 1000*resdiff{gather}./pick_times;
      resdiff_signed_norm{gather} = 1000*resdiff_signed{gather}./pick_times;
      
      % exclude picks with large resid error from 200 samples down
      
%       ttimes2 = ttimes(resdiff<=median(resdiff),:);
%       ini_tv2 = ini_tv(resdiff<=median(resdiff),:);
%       scale_fn2 = scale_fn(resdiff<=median(resdiff),:);

      v_est2{gather} = v_est{gather}((resdiffnorm{gather}<=0.5) | (t_est{gather}<400));
      t_est2{gather} = t_est{gather}((resdiffnorm{gather}<=0.5) | (t_est{gather}<400));

      velout(:,gather,2) = interp1([0;t_est2{gather}],[v_est2{gather}(1);v_est2{gather}],2:2:3000);
      
      
%       [v_est{gather} t_est{gather}] = invert_T_V_from_picks(ttimes2,offsets(:,gather),ini_tv2);
%       
%        t_est{gather}(1)=t_est{gather}(1)+1;
%       velout(:,gather,3) = interp1([0;t_est{gather}],[v_est{gather}(1);v_est{gather}],2:2:3000);
     
      
%   end
  
  
  tvpairs{gather}(:,1) = t_est{gather};
  tvpairs{gather}(:,2) = v_est{gather};
  
  tvpairs2{gather}(:,1) = t_est2{gather};
  tvpairs2{gather}(:,2) = v_est2{gather};

end

% write some output


stack_traces{2,2} = stack;
stack_traces{4,2} = velout;
i_block = str2double(i_block);
node_segy_write(stack_traces,i_block,seismic.s_rate/1000,outdir);

save(strcat(outdir,'rms_picks_block',int2str(i_block),'.mat'),'ntraces','tvpairs','tvpairs2',...
    'stack_traces','stack','stack_mask','resdiff','resdiffnorm','resdiff_signed','resdiff_signed_norm');


                        end

function vint = velfun(vrms,ttimes)
tts = [0;ttimes(1:end-1)];
vrms_s = [vrms(1);vrms(1:end-1)];

vint = sqrt(((vrms.^2.*ttimes) - (vrms_s.^2.*tts)) ./ (ttimes-tts));
end

