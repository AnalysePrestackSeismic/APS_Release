function [] = register_seis_vols(meta_path,iblock,output_dir,start_time,end_time,reg_type,reg_smth)
%
% register_vols: register two seismic volumes
%
% The input volume should be two fold gathers with the base volume in the
% first trace and the monitor in the second
%
% Arguments:
%   meta_path:      pathname of meta file for 2-fold gathers
%   output_dir:     path to output directory. Must end with "/".
%   start_time:     start time in ms
%   end_time:       end_time in ms       
%   reg_type:       'affine' or 'nonrigid' registration (defaults to nonrigid)
%   reg_smth:       fluid sigma factor, or smoothness (default is 4)
%   iblock:         block number to process
%
% Writes to disk
%   base_vol_reg:   registered base volume
%   regx:           forward (base to monitor) transformation x-component
%   regy:           forward (base to monitor) transformation y-component
%   regz:           forward (base to monitor) transformation z-component
%
% segy_block_size should be called beforehand to set the block size and
% overlap for the analysis
%
% The "register_volumes" function threads internally so only one instance
% can be run per node.

reg_smth = str2double(reg_smth);
start_time = str2double(start_time);
end_time = str2double(end_time);

job_meta = load(meta_path);

start_samp = 1+1000*start_time/job_meta.s_rate;
end_samp = 1+1000*end_time/job_meta.s_rate;
n_samps = 1+end_samp-start_samp;


[block_meta, traces, ilxl, offsets] = node_segy_read(meta_path,'1',iblock);

min_pkey = min(ilxl(:,1));
max_pkey = max(ilxl(:,1));
min_skey = min(ilxl(:,2));
max_skey = max(ilxl(:,2));

n_pkeys = 1+max_pkey-min_pkey/job_meta.pkey_inc;
n_skeys = 1+max_skey-min_skey/job_meta.skey_inc;

% initialise the volumes to load the traces into

base_vol = zeros(n_samps,n_skeys*n_pkeys);
monitor_vol = base_vol;

base_traces=(traces(start_samp:end_samp,1:2:end-1));
monitor_traces=(traces(start_samp:end_samp,2:2:end));

trace_index = 1+(ilxl(1:2:end-1,2)-min_skey)./job_meta.skey_inc + ((ilxl(1:2:end-1,1)-min_pkey)./job_meta.pkey_inc).*n_skeys;

base_vol(:,trace_index) = base_traces;
monitor_vol(:,trace_index) = monitor_traces;

base_vol = double(reshape(base_vol,n_samps,n_skeys,[])); % convert to double for "register volumes" function
monitor_vol = double(reshape(monitor_vol,n_samps,n_skeys,[])); % convert to double for "register volumes" function


if strcmpi(reg_type,'affine')
    
    % affine registration
    
    [base_vol_reg,regz,regx,regy] = register_volumes(base_vol,monitor_vol,struct('Interpolation','Cubic','Registration','Affine'));
    
else
    
    % non-rigid registration
    
    [base_vol_reg,regz,regx,regy] = register_volumes(base_vol,monitor_vol,struct('Interpolation','Cubic','Registration','NonRigid','SigmaFluid',reg_smth));
    
end

% for the output
segy_out{1,1} = 'Meta data for output files';
segy_out{1,2}{1,1} = ilxl(2:2:end,:);
segy_out{1,2}{2,1} = int32(zeros(block_meta.n_traces,1));
segy_out{1,1} = sprintf('%-3200.3200s','Image registration');
segy_out{1,3} = 'is_gather'; % 1 is yes, 0 is no

% output_dir = '/data/EGY/segy/WDDM/2014/Working/registration/';
% check to see if the directory exists
if exist(output_dir,'dir') == 0
    mkdir(output_dir);
end

reg_smth = num2str(reg_smth);

% only output traces for which there was input:

segy_out{2,1} = strcat('base_reg_',reg_type,'_smth',reg_smth);
segy_out{2,2} = base_vol_reg(:,trace_index);
segy_out{2,3} = 0;

segy_out{3,1} = strcat('regx_',reg_type,'_smth',reg_smth);
segy_out{3,2} = regx(:,trace_index);
segy_out{3,3} = 0;

segy_out{4,1} = strcat('regy_',reg_type,'_smth',reg_smth);
segy_out{4,2} = regy(:,trace_index);
segy_out{4,3} = 0;

segy_out{5,1} = strcat('regz_',reg_type,'_smth',reg_smth);
segy_out{5,2} = regz(:,trace_index);
segy_out{5,3} = 0;

node_segy_write(segy_out,str2double(iblock), block_meta.s_rate, output_dir);

end
