i_body = 1;
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
area_thresh = 10000;
pkeys = job_meta.pkey_min:job_meta.pkey_inc:job_meta.pkey_max;
skeys = job_meta.skey_min:job_meta.skey_inc:job_meta.skey_max;

for ii = 1:1:CC.NumObjects
%    I_vol(CC.PixelIdxList{ii}) = STATS(ii).Area; % volume
    
    if STATS(ii).Area > area_thresh
        z_col = 2;
        skey_col = 1;
        pkey_col = 3;
      
        pos = STATS(ii).PixelList;
        sort_z_as = sortrows(pos,z_col);
        sort_z_ds = sortrows(pos,-z_col);
        pos_uniq = unique(pos(:,[pkey_col skey_col]),'rows');
        % col 1 skey
        % col 2 z
        % col 3 pekyn           
        pos_body_surface_top(:,1) = pkeys(pos_uniq(:,1));
        pos_body_surface_top(:,2) = skeys(pos_uniq(:,2));
        
        [~,Locb] = ismember(pos_uniq,sort_z_ds(:,[pkey_col skey_col]),'rows');
        pos_body_surface_top(:,3) = (sort_z_ds(Locb,2)-1).*job_meta.s_rate/1000;
        
        pos_body_surface_bot(:,1) = pkeys(pos_uniq(:,1));           
        pos_body_surface_bot(:,2) = skeys(pos_uniq(:,2));        
        [~,Locb] = ismember(pos_uniq,sort_z_as(:,[3 1]),'rows');
        pos_body_surface_bot(:,3) = (sort_z_as(Locb,2)-1).*job_meta.s_rate/1000;
        
        top_asc = strcat(job_meta.output_dir,'/surface_body_',num2str(i_body),'_top.txt');
        bot_asc = strcat(job_meta.output_dir,'/surface_body_',num2str(i_body),'_bot.txt');
        
        dlmwrite(top_asc,pos_body_surface_top,'delimiter','\t','precision', 8)
        dlmwrite(bot_asc,pos_body_surface_bot,'delimiter','\t','precision', 8)
        
        i_body = i_body + 1;
        
        clear pos sort_z_as sort_z_ds pos_uniq Locb pos_body_surface_top pos_body_surface_bot
    end
end