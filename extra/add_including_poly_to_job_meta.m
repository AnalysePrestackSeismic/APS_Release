function [] = add_including_poly_to_job_meta(job_meta_path,ils,xls)
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
%   ils: Array of inline numbers
%   xls: Array of xline numbers
% Saves to Disk:
%   new job meta file with cropped live blocks
% Note: The blocks that are fully inside the polygon will qualify
job_meta.ply.il=ils;
job_meta.poly.xl=xls;
poly.method=1;
ils=double(ils);
xls=double(xls);
ils=ils';
xls=xls';


inline_sort_blocks(job_meta_path)
job_meta = load(job_meta_path);                 % Load job meta information
counter=1;

for i_block=1:length(job_meta.liveblocks);
    blk=job_meta.liveblocks(i_block);
    ilxl_blk=job_meta.block_keys(blk,:);
    
    % Define the corner points of the block
    il_corner_pts= [ilxl_blk(1) ilxl_blk(2)];
    xl_corner_pts= [ilxl_blk(3) ilxl_blk(4)];
    
    % Check if the corner points are inside deined polygon
    [IN ON]= inpolygon(il_corner_pts,xl_corner_pts,ils,xls);
    
    %If both the corner points are inside the polygon this block is
    %considered as a likve block
    
    if sum(IN)==2
       liveblocks_poly(counter)=job_meta.liveblocks(i_block);
       counter=counter+1;
    end
end
job_meta.liveblocks = liveblocks_poly';
str_date = date;
str_date = regexprep(str_date, '-', '');
job_meta_path_poly=strcat(job_meta.output_dir,'job_meta/','job_meta_poly_',str_date,'.mat');

save(job_meta_path_poly,'-struct','job_meta','-v7.3');
end
