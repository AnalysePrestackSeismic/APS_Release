function [ ] = wrapper_digi_on_inline( iL,job_meta_path,startvol,volinc,endvol,maxzout,wavevar)
% Wrapper function for running digi on one inline

inline = str2double (iL);
job_meta = load(job_meta_path);

flag_block_keys = find(job_meta.block_keys(job_meta.liveblocks,1)<inline & job_meta.block_keys(job_meta.liveblocks,2)>inline); % FInd the blocks

for blk = 1: size(flag_block_keys)
    iblock = num2str (flag_block_keys(blk));
   int_grad_inv_proj(job_meta_path,iblock,startvol,volinc,endvol,iL,maxzout,wavevar)
end

end
