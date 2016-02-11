function [] = crop_live_blocks(job_meta_path,il1,xl1,il2,xl2)
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
%% Function to crop thelive blocks entry in job meta file by a rectangle defined by inlines and crosspline rages
% input: 
%   job_meta_path : path to job meta file
%   il1 : smallest inline
%   xl1 : smallest cross-line
%   il2 : largest inline
%   xl2 :largest cross line

% output: void

% writes to disk: New job meta file with cropped live blocks

% Displays: The patn of new job meta file

%%
inline_sort_blocks(job_meta_path)
job_meta = load(job_meta_path);                 % Load job meta information
counter=1;

% Loop through all the live blocks in job meta files
for i_block=1:length(job_meta.liveblocks);
    blk=job_meta.liveblocks(i_block);
    ilxl_blk=job_meta.block_keys(blk,:);
    % check if the live block is with in the user defined boundaries
    if (ilxl_blk(1)>il1 && ilxl_blk(2)<il2 && ilxl_blk(3)>xl1 && ilxl_blk(4)<xl2)
       liveblocks_crop(counter)=job_meta.liveblocks(i_block);
       counter=counter+1;
    end
end
job_meta.liveblocks = liveblocks_crop';
str_date = date;
str_date = regexprep(str_date, '-', '');
job_meta_path_crop=strcat(job_meta.output_dir,'job_meta/','job_meta_crop',str_date,'.mat');

save(job_meta_path_crop,'-struct','job_meta','-v7.3');
disp(job_meta_path_crop);
end
