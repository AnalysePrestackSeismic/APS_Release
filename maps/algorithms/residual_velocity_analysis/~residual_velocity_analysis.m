function [] = residual_velocity_analysis( gather_meta_path,velocity_meta_path,i_block )

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
stack_otm = 30;
pick_otm = 40;


% ============================================================================================================
%
% load the data
%

gather_meta = load(gather_meta_path);

velocity_meta = load(velocity_meta);

% read the gathers and velocities in from segy

% trace data is ordered by cdp, offset, sample

% incoming trace data is one stream of traces for efficiency
% so reshape it into a 3 dimensional array for processing

[seismic, traces{3,2}, traces{1,2}{1,1}, offsets] = node_segy_read(gather_meta_path,'1',i_block);
seismic.fold = max(seismic.trace_ilxl_bytes(:,7));
traces{3,2} = reshape(traces{3,2},size(traces{3,2},1),seismic.fold,[]);

[velocity, vel_traces{3,2}, vel_traces{1,2}{1,1},~] = node_segy_read(velocity_meta_path,'0',i_block)

seismic.n_gathers = size(traces{3,2},3);
    
% check to make sure it read something if not exit
if isempty(traces{3,2}) == 1 &&  isempty(traces{1,2}{1,1}) == 1 && isempty(offsets) == 1
    return
end
if isempty(vel_traces{3,2}) == 1 &&  isempty(vel_traces{1,2}{1,1}) == 1 
    return
end

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
traces{2,1} = strcat(filename,'_trim_shifts_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
traces{3,1} = strcat(filename,'_nmo_data_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
traces{2,3} = 1;
traces{3,3} = 1;
traces{2,2} = zeros(in_n_samples,n_traces_gather,n_traces,'single');


% stack result
filename2 = 'stack';
stack_traces{2,1} = strcat(filename2,'_input_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
stack_traces{2,2} = zeros(in_n_samples,n_traces,'single');
stack_traces{3,1} = strcat(filename2,'_output_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
stack_traces{3,2} = zeros(in_n_samples,n_traces,'single');
%stack_traces{4,1} = strcat(filename2,'_posttrim_no_resid_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
stack_traces{4,1} = strcat(filename2,'_velocity_',num2str(startvol),'-',num2str(endvol),'_',num2str(zsmooth));
stack_traces{4,2} = zeros(in_n_samples,n_traces,'single');
stack_traces{2,3} = 0;
stack_traces{3,3} = 0;
stack_traces{4,3} = 0;

% ============================================================================================================
% calculate offsets from angles and velocities and make mute masks

angles = [stack_itm stack_otm pick_otm];

stack_mute=zeros(size(vel_traces{3,2});
pick_mute=zeros(size(vel_traces{3,2});
stack_fold=zeros(velocity.nsamples,velocity.n_gathers);
pick_fold=zeros(velocity.nsamples,velocity.n_gathers);

for gather_idx=1:seismic.ngathers
    
    offset_lookup = offset_vs_angle(velocity_meta,angles,vel_traces{3,2});
    
    for off_idx=1:seismic.fold
        stack_mute(:,off_idx,gather_idx) = (offsets(gather_idx) >= offset_lookup(:,1,gather_idx)) ...
            .* (offsets(gather_idx) <= offset_lookup(:,2,gather_idx));
        pick_mute(:,off_idx,gather_idx) = (offsets(gather_idx) <= offset_lookup(:,3,gather_idx));
    end
    stack_fold(:,gather_idx) = sum(stack_mute,2);
    pick_fold(:,gather_idx) = sum(pick_mute,2);
end

% pad and smooth mutes then remove padding

smth = [1;1;1;1;1];
stack_mute = [stack_mute;repmat(stack_mute(end),size(smth,1))];
stack_mute = conv2(stack_mute,smth,'same')./sum(smth);
stack_mute = stack_mute(:velocity.nsamples,:,:);


% ============================================================================================================
% make a stack for picking freq content and waterbottom

% interp mute mask onto seismic sampling

ss = seismic.srate;
sn = seismic.nsamples;
vs = velocity.srate;
vn = velocity.nsamples;

% pad mute with extra row at start and interpolate:
stack_mute_interp = interp1([ss vs:vs:vs*vn],[stack_mute(1,:,:);stack_mute];,ss:ss:ss*sn);

clear ss; clear sn; clear vs; clear vn;

stack_traces{2,2} = (sum(traces{3,2},2);


end

