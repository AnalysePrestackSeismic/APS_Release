function [] = copy_job_meta(copy_from,copy_to,varargin)
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

