function [] = add_mute_to_job_meta(job_meta_path,method,varargin)
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
%% Funtion to add mute to job meta file. When you want to mute the data before doing mica (especially useful for offset gathers)

%   INPUTS:
%     Method = '1' if u want give mute as a set (two arrays ) of time and offset
%     Method = '2' if you want to give mute as text file
%   OUPUTS: 
% WRITES TO DISK:
%   saves wb_path in job meta file

%%
job_meta = load(job_meta_path);   % load job meta file

method=str2double(method);


if method==1
    sample_number = varargin{1}/(job_meta.s_rate/1000);
    sample_number = [sample_number job_meta.n_samples{1,1}];
    offset = varargin{2};
    offset= [offset offset(end)];
    
    job_meta.mute=[sample_number' offset'];
end    

 
save(job_meta_path,'-struct','job_meta','-v7.3');  
  
end
