function [] = add_mute_to_job_meta(job_meta_path,method,varargin)
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
