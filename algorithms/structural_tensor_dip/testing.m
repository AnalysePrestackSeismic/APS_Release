% Test Structural Tensor Dip
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
job_meta_path = '/data/TZA/dtect/2014_tachuii_presdm_kirchoff/Misc/output_matlab/job_meta/job_meta_20Mar2015.mat';

i_block = '255';
vol_index = '1';
sigma = '3';
scale_sigma = '4';
start_slab = '1';
end_slab = '2000';

% set output path
aperture = num2str(str2num(scale_sigma)*str2num(sigma));
job_meta = load(job_meta_path);
save([job_meta_path, '_no_aperture'],'-struct','job_meta','-v7.3');
add_aperture_to_job_meta(job_meta_path,aperture);

structural_tensor_dip(job_meta_path,i_block,vol_index,sigma,scale_sigma,start_slab,end_slab)

% Remove Aperture
add_aperture_to_job_meta(job_meta_path,['-',aperture]);