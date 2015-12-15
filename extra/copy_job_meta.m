function [] = copy_job_meta(copy_from,copy_to,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

  job_meta_from = load(copy_from);
  job_meta_to = load(copy_to);
  
  fields_from = fieldnames(job_meta_from);
  
  for i_field = 1:size(fields_from,1)
        fprintf('Checking field, %s\n',fields_from{i_field});
        if ~isfield(job_meta_to, fields_from{i_field})
            job_meta_to.(fields_from{i_field}) = job_meta_from.(fields_from{i_field});
        end
  end
  
  save(copy_to,'-struct','job_meta_to','-v7.3');  

  
end

