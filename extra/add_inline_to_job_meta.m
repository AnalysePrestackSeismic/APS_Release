function [] = add_inline_to_job_meta(job_meta_path,il)
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
%% Function to add a polygon to job meta file and crop the live blocks by the polygon
% INPUT:
%   job_meta_path : path of original job meta file
%   il : Inline Number
% Saves to Disk:
%   new job meta file with live blocks that have this inline
%%
il=str2double(il);

inline_sort_blocks(job_meta_path)
job_meta = load(job_meta_path);                                     % Load job meta information

x=(job_meta.block_keys(:,2)>=il).*(job_meta.block_keys(:,1)<il);    % find the blocks with this inline
ind=find(x);
new_liveblocks=ind(ismember(ind,job_meta.liveblocks));              % check if the selcted blocks are live blocks or not
job_meta.liveblocks=new_liveblocks;                                 % update live blocks in job meta file

str_date = date;
str_date = regexprep(str_date, '-', '');
job_meta_path_poly=strcat(job_meta.output_dir,'job_meta/','job_meta_il_',num2str(il),'_',str_date,'.mat');

save(job_meta_path_poly,'-struct','job_meta','-v7.3');              %save new job meta file
end
