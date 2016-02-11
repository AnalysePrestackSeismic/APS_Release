function [ ] = wrapper_digi_on_inline( iL,job_meta_path,startvol,volinc,endvol,maxzout,wavevar)
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
% Wrapper function for running digi on one inline

inline = str2double (iL);
job_meta = load(job_meta_path);

flag_block_keys = find(job_meta.block_keys(job_meta.liveblocks,1)<inline & job_meta.block_keys(job_meta.liveblocks,2)>inline); % FInd the blocks

for blk = 1: size(flag_block_keys)
    iblock = num2str (flag_block_keys(blk));
   int_grad_inv_proj(job_meta_path,iblock,startvol,volinc,endvol,iL,maxzout,wavevar)
end

end
