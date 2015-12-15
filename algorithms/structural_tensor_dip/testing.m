% Test Structural Tensor Dip

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