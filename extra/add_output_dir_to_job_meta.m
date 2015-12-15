function [] = add_output_dir_to_job_meta(job_meta_path,output_path)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

  job_meta = load(job_meta_path);
   
  job_meta.output_dir = output_path;
  
  job_meta_path = strsplit(job_meta_path,'.mat');
  job_meta_path = cell2mat(job_meta_path);
  job_meta_path = [job_meta_path,'_outputdir.mat'];    
  save(job_meta_path,'-struct','job_meta','-v7.3');      

end

