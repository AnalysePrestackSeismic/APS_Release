function [] = residual_velocity_analysis_aruba_test( gather_meta_path,i_block,outdir,apply_nmo_switch)
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
% apply_nmo_switch:     should be set to 1 if input gathers do not have nmo applied           
%
% velocity field should be sampled on same x,y grid as gathers but can be
% decimated in t
%

% to reduce printout in compilied version turned all warning off
warning off all;
apply_nmo_switch = str2double(apply_nmo_switch);

% ============================================================================================================
% hardwired parameters

stack_itm = 0;
stack_otm = 20;
pick_otm = 20;

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

% manual picks we are attempting to replicate:

%       [[875.1       , 1508        ],
%        [966.6       , 1517        ],
%        [1153        , 1571        ],
%        [1323        , 1615        ],
%        [1525        , 1651        ],
%        [1594        , 1633        ],
%        [1740        , 1687        ],
%        [1878        , 1687        ],
%        [2106        , 1804        ],
%        [2314        , 1840        ],
%        [2459        , 1894        ]]
% 

[velstr velilxl veltrace] =  segy_to_mat('189','193','/data/ABW/segy/3D_2014_PSDM/BG_internal/velocity_test/RMS_Velocities_block846.segy');

veltimes = 0.001*((velstr.s_rate/1000):(velstr.s_rate/1000):(velstr.s_rate/1000)*velstr.n_samples)';



% picktimes = 0.001*[875;976;1153;1323;1525;1594;1740;1878;2106;2314;2459];

% veltrace = interp1([0 6000],[1500 1500],100:100:5000);
% 
% veltimes = ([100:100:5000]')./1000;
% 
% veltrace = veltrace';


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

if isempty(traces{3,2}) == 1 &&  isempty(traces{1,2}{1,1}) == 1 && isempty(offsets) == 1
    return
end
%if isempty(vel_traces{3,2}) == 1 &&  isempty(vel_traces{1,2}{1,1}) == 1 
%    return
%end

% veltrace = repmat(veltrace,1,seismic.n_gathers);
% veltimes = repmat(veltimes,1,seismic.n_gathers);

gather_ilxl = unique(traces{1,2}{1,1},'rows');

min_il = double(traces{1,2}{1,1}(1,1));
max_il = double(traces{1,2}{1,1}(end,1));

min_xl = double(traces{1,2}{1,1}(1,2));
max_xl = double(traces{1,2}{1,1}(end,2));

il_inc = double(gather_meta.pkey_inc);
xl_inc = double(gather_meta.skey_inc);

num_xls = 1 + (max_xl-min_xl)/xl_inc;
num_horz = 11;

picktimes = zeros(seismic.n_gathers,num_horz);

for horz = 1:num_horz
    horz_file = strcat('/data/ABW/dtect/Aruba_DanB_Time/Misc/veltest_h',num2str(horz),'.txt');
    horz_in = dlmread(horz_file);
    il_idx = int32((horz_in(:,1)-min_il)./il_inc);
    xl_idx = int32((horz_in(:,2)-min_xl)./xl_inc);
    ilxl_idx = il_idx.*num_xls + xl_idx + 1;
    picktimes(ilxl_idx,horz) = horz_in(:,3);
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
% tt = ((1:nsamples).*100)';

%tts = [0;tt(1:end-1)];
%veltrace_s = [veltrace(1);veltrace(1:end-1)];

%vint = sqrt(((veltrace.^2.*tt) - (veltrace_s.^2.*tts)) ./ (tt-tts));

% this is how to do it if you have more than one vel trace but only one
% column of times

vint = bsxfun(@velfun,veltrace,veltimes);



% vint = sqrt(bsxfun(@rdivide,(bsxfun(@times,(veltrace.^2),tt) - bsxfun(@times,(veltrace_s.^2),tts)),(tt-tts)));



% =====================================================================

% angle calculation

% sine theta = (offset * vint) / (t * vrms^2)
% => offset = (t * vrms^2 * sine theta) / vint

% tan theta = 0.5 * offset / depth = 0.5 * offset / (0.5 * twtt * vrms)
% => offset = tan theta * twtt * vrms

for ang_idx=1:max(size(angles))
    
    offset_traces(:,ang_idx,:) = bsxfun(@times,veltimes,veltrace * tan(angles_rad(ang_idx)));
    
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

vs=velstr.s_rate/1000;
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

%vel_out = zeros(size(picktimes));
times_out = zeros(size(picktimes));

% we want output to be inline,xline,time,velocity
gap = 12;
output_matrix = zeros((seismic.n_samples/gap)*seismic.n_gathers,4);
outcount = 0;
gathcount = 0;
% run just for a subset for testing
for gather = 1:seismic.n_gathers
    
%    if sum(picktimes(gather,:))>0
        
        
        % this code uses horizon times from picktimes matrix ==========
%         ptimes = sort(picktimes(gather,picktimes(gather,:)>0))';
%         pick_idx = int32(1000*ptimes./gather_meta.s_rate);
        % ==============================================================

        [trim_shifts,~,~,trim_coeffs] = trim_shifts_calculate('calc',gather_meta,traces{3,2}(:,:,gather),'1','4');
        
        % alternatively use xcor coefficent peaks =====================
        %
        coeff_peaks_out = zeros(size(trim_coeffs));
        while max(trim_coeffs)>0.8
            
            [~,peak_idx] = max(trim_coeffs);
            coeff_peaks_out(peak_idx) = trim_coeffs(peak_idx);
            start_mask=max(peak_idx-gap,1);
            end_mask=min(peak_idx+gap,seismic.n_samples);
            
            trim_coeffs(start_mask:end_mask) = 0;
        end
        pick_idx = find(coeff_peaks_out);
        ptimes = pick_idx * seismic.s_rate/1000;

        %==============================================================
        
        numpicks = size(ptimes,1);
        outcount = outcount + numpicks;
        vels = 0.001*interp1(veltimes*1000,veltrace(:,gather),ptimes); % times are in ms, so make vels m/ms
%         times_out(gather,1:numpicks)=ptimes;

        t0times = bsxfun(@minus,trim_shifts(pick_idx,:),ptimes);
        %t0times2 = bsxfun(@plus,-trim_shifts(pick_idx,:,gather),pick_times);
        %t0times2 = bsxfun(@minus,trim_shifts(pick_idx,:,gather),pick_times);
        ttimes = sqrt(t0times.^2 + bsxfun(@rdivide,double(offsets(:,gather)'.^2),vels.^2));
        
        % set times in mute zone to be NaN
        
        ttimes(stack_mute(pick_idx,:,gather)<1) = NaN;
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
    
    gathcount = gathcount + 1;
    
    if mod(gathcount,10) == 0
        disp(['Done gather ',num2str(gather)]);
    end
    
%    end
        
end

dlmwrite(strcat(outdir,'veltest_rms_picks_xcor_peaks_gap12_mute20_block846.txt'),output_matrix(1:outcount,:));


  

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

