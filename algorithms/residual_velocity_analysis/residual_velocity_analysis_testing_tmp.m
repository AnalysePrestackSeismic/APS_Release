function [] = residual_velocity_analysis_testing( gather_meta_path,velocity_meta_path,i_block )

% Automatic residual velocity analysis
%   
% Inverts residual time errors on NMO corrected gathers to give updated
% velocity field
%
% gather_meta_path:     pathname to .mat file for the NMO-corrected gathers
% velocity_meta_path:   pathname to .mat file for the initial velocities
% i_block:              block number to analyse
%
% velocity field should be sampled on same x,y grid as gathers but can be
% decimated in t
%

% to reduce printout in compilied version turned all warning off
warning off all;
% ============================================================================================================
% hardwired parameters

stack_itm = 5;
stack_otm = 60;
pick_otm = 70;


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
veltrace = repmat(veltrace,1,seismic.n_gathers);
% check to make sure it read something if not exit
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

offsets = reshape(offsets,seismic.fold,[]);

% compare trace offsets with those calculated for the mute angles and
% create a mute mask for each gather

mute1 = squeeze(offset_traces(:,1,:));
mute2 = squeeze(offset_traces(:,2,:));
mute3 = squeeze(offset_traces(:,3,:));

% make matrix to keep n traces live when applying otm

otm_limit = [ones(4,1);zeros(seismic.fold-4,1)];
   

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

traces{3,2} = traces{3,2}.*stack_mute;

% run just for a subset for testing

[trim_shifts,~,~,trim_coeffs] = trim_shifts_calculate('calc',gather_meta,traces{3,2}(:,:,1:100),'1','1');

% pick xcorrelation peaks - these are times to pick velocities at

% fairly sure this is v inefficient way of doing it
% also try smoothing the coefficents to get more consistent time picks

coeff_peaks_out = zeros(size(trim_coeffs));

smth = [1 2 3 4 5 6 7 8 9 8 7 6 5 4 3 2 1];
smth = smth / sum(smth);
pad = ceil(size(smth,2)/2);

coeffs_smth = convn([repmat(trim_coeffs(:,1),1,pad) trim_coeffs repmat(trim_coeffs(:,size(trim_coeffs,2)),1,pad)],smth,'same');
coeffs_smth = coeffs_smth(:,(pad+1):(size(coeffs_smth,2)-pad));

for trc = 1:size(trim_coeffs,2);
    
    while max(trim_coeffs(:,trc))>0.8
        
        [~,peak_idx] = max(trim_coeffs(:,trc));
        coeff_peaks_out(peak_idx,trc) = trim_coeffs(peak_idx,trc);
        start_mask=max(peak_idx-50,1);
        end_mask=min(peak_idx+50,seismic.n_samples);
        
        trim_coeffs(start_mask:end_mask,trc) = 0;
        
    end
    
end

for trc = 1:size(coeffs_smth,2);
    
    while max(coeffs_smth(:,trc))>0.8
        
        [~,peak_idx] = max(coeffs_smth(:,trc));
        coeff_sm_peaks_out(peak_idx,trc) = coeffs_smth(peak_idx,trc);
        start_mask=max(peak_idx-50,1);
        end_mask=min(peak_idx+50,seismic.n_samples);
        
        coeffs_smth(start_mask:end_mask,trc) = 0;
        
    end
    
end


% ============================================================================================================

% calc travel times for trim picks (assume straight ray)

% tt = sqrt(t0^2 + (offset^2/v^2))

%for gather = 1:seismic.n_gathers

velout = zeros(1500,100);
velout_sm = zeros(1500,100);

for gather = 1:100
    
  pick_idx = find(coeff_peaks_out(:,gather));
  pick_times = pick_idx * seismic.s_rate/1000;
  
  vels = 0.001*interp1([0 100:100:5000],[veltrace(1,gather);veltrace(:,gather)]',pick_times);
  
  t0times = bsxfun(@minus,trim_shifts(pick_idx,:,gather),pick_times);
  %t0times2 = bsxfun(@plus,-trim_shifts(pick_idx,:,gather),pick_times);
  %t0times2 = bsxfun(@minus,trim_shifts(pick_idx,:,gather),pick_times);
  ttimes = sqrt(t0times.^2 + bsxfun(@rdivide,double(offsets(:,gather)'.^2),vels.^2));

  % set times in mute zone to be NaN
  
  ttimes(stack_mute(pick_idx,:,1)<1) = NaN;
  
  [v_est{gather} t_est{gather}] = invert_T_V_from_picks(ttimes,offsets(:,gather));
  
  % hardwire first vel to 1540
  
  % v_est{gather}(1)=1.54;
  t_est{i}(1)=t_est{i}(1)+1;
  velout(:,gather) = interp1([0;t_est{gather}],[1.54;v_est{gather}],2:2:3000);
end

for gather = 1:100
    
  pick_idx = find(coeff_sm_peaks_out(:,gather));
  pick_times = pick_idx * seismic.s_rate/1000;
  
  vels = 0.001*interp1([0 100:100:5000],[veltrace(1,gather);veltrace(:,gather)]',pick_times);
  
  t0times = bsxfun(@minus,trim_shifts(pick_idx,:,gather),pick_times);
  %t0times2 = bsxfun(@plus,-trim_shifts(pick_idx,:,gather),pick_times);
  %t0times2 = bsxfun(@minus,trim_shifts(pick_idx,:,gather),pick_times);
  ttimes = sqrt(t0times.^2 + bsxfun(@rdivide,double(offsets(:,gather)'.^2),vels.^2));

  % set times in mute zone to be NaN
  
  ttimes(stack_mute(pick_idx,:,1)<1) = NaN;
  
  [v_sm_est{gather} t_sm_est{gather}] = invert_T_V_from_picks(ttimes,offsets(:,gather));
  
  % hardwire first vel to 1540
  
  % v_est{gather}(1)=1.54;
  
  velout_sm(:,gather) = interp1([0;t_sm_est{gather}],[1.54;v_sm_est{gather}],2:2:3000);
  

end



end

function vint = velfun(vrms,ttimes)
tts = [0;ttimes(1:end-1)];
vrms_s = [vrms(1);vrms(1:end-1)];

vint = sqrt(((vrms.^2.*ttimes) - (vrms_s.^2.*tts)) ./ (ttimes-tts));
end

