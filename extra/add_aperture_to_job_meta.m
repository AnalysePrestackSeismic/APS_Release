function [job_meta_path] = add_aperture_to_job_meta(job_meta_path,aperture)
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

