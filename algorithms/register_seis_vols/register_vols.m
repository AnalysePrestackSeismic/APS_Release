function [base_vol_reg regx regy regz] = register_vols(meta_path,start_time,end_time,reg_type,reg_smth,iblock)
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
%
% register_vols: register two seismic volumes
%
% The input volume should be two fold gathers with the base volume in the
% first trace and the monitor in the second
%
% Arguments:
%   meta_path:      pathname of meta file for 2-fold gathers
%   start_time:     start time in ms
%   end_time:       end_time in ms       
%   reg_type:       'affine' or 'nonrigid' registration (defaults to nonrigid)
%   reg_smth:       fluid sigma factor, or smoothness (default is 4)
%   iblock:         block number to process
%
% Outputs:
%   base_vol_reg:   registered base volume
%   regx:           forward (base to monitor) transformation x-component
%   regy:           forward (base to monitor) transformation y-component
%   regz:           forward (base to monitor) transformation z-component
%
% segy_block_size should be called beforehand to set the block size and
% overlap for the analysis
%


job_meta = load(meta_path);

start_samp = 1+1000*start_time/job_meta.s_rate;
end_samp = 1+1000*end_time/job_meta.s_rate;
n_samps = 1+end_samp-start_samp;


[block_meta, traces, ilxl, offsets] = node_segy_read(meta_path,'1',num2str(iblock));

n_skeys = 1+max(ilxl(:,2))-min(ilxl(:,2))/job_meta.skey_inc;

base_vol = double(reshape(traces(start_samp:end_samp,1:2:end-1),n_samps,n_skeys,[]));
monitor_vol = double(reshape(traces(start_samp:end_samp,2:2:end),n_samps,n_skeys,[]));

if strcomp(reg_type,'affine')

% affine registration

[base_vol_reg,regz,regx,regy] = register_volumes(base_vol,monitor_vol,struct('Interpolation','Cubic','Registration','Affine'));

else

% non-rigid registration 

[base_vol_reg,regz,regx,regy] = register_volumes(base_vol,monitor_vol,struct('Interpolation','Cubic','Registration','NonRigid','SigmaFluid',reg_smth));

end

end
