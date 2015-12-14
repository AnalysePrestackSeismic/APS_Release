function [] = add_aperture_to_job_meta(job_meta_path,aperture)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

  job_meta = load(job_meta_path);
   
  job_meta.output_dir = output_path;
    
  save(job_meta_path,'-struct','job_meta','-v7.3');  
  
end

