function [] = remove_aperture_to_job_meta(job_meta_path)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

  job_meta = load(job_meta_path);
  
  job_meta.block_keys(:,1) = job_meta.block_keys(:,1)+job_meta.aperture;
  job_meta.block_keys(:,2) = job_meta.block_keys(:,2)-job_meta.aperture;
  job_meta.block_keys(:,3) = job_meta.block_keys(:,3)+job_meta.aperture;
  job_meta.block_keys(:,4) = job_meta.block_keys(:,4)-job_meta.aperture;
  
  job_meta = rmfield(job_meta,'aperture');
  save(job_meta_path,'-struct','job_meta','-v7.3');  
  
end

