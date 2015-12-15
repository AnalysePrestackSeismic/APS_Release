function [job_meta_path] = add_aperture_to_job_meta(job_meta_path,aperture)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

  job_meta = load(job_meta_path);
  
  job_meta.block_keys(:,1) = job_meta.block_keys(:,1)-str2num(aperture);
  job_meta.block_keys(:,2) = job_meta.block_keys(:,2)+str2num(aperture);
  job_meta.block_keys(:,3) = job_meta.block_keys(:,3)-str2num(aperture);
  job_meta.block_keys(:,4) = job_meta.block_keys(:,4)+str2num(aperture);
   
  job_meta.aperture = str2num(aperture);
  job_meta_path = strsplit(job_meta_path,'.mat');
  job_meta_path = cell2mat(job_meta_path);
  job_meta_path = [job_meta_path,'_aperture.mat'];    
  save(job_meta_path,'-struct','job_meta','-v7.3');  
  
end

