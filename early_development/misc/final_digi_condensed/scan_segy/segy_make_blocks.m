function [block_keys,n_blocks] = segy_make_blocks(job_meta_path,pkey_blocks,skey_blocks,tkey_blocks)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% if gathers

    job_meta = load(job_meta_path);

    pkey_blocks = str2num(pkey_blocks);
    skey_blocks = str2num(skey_blocks);
    tkey_blocks = str2num(tkey_blocks);

    % Need to correct
    pkey_inc = floor((mode(job_meta.pkey_max) - mode(job_meta.pkey_min))/pkey_blocks);
    skey_inc = floor((mode(job_meta.skey_max) - mode(job_meta.skey_min))/skey_blocks);

    pkeys_min = job_meta.pkey_min:pkey_inc:job_meta.pkey_max-pkey_inc;
    skeys_min = job_meta.skey_min:skey_inc:job_meta.skey_max-skey_inc;
    pkey_blocks = size(pkeys_min,2);
    skey_blocks = size(skeys_min,2);
    
    pkeys_max = pkeys_min - 1;
    pkeys_max = pkeys_max(2:end);
    pkeys_max(end+1) = mode(job_meta.pkey_max);

    skeys_max = skeys_min - 1;
    skeys_max = skeys_max(:,2:end);  
    skeys_min = repmat(skeys_min,pkey_blocks,1);
    skeys_max = repmat(skeys_max,pkey_blocks,1);
    skeys_max(:,end+1) = repmat(mode(job_meta.skey_max),pkey_blocks,1);
    pkeys_min = repmat(pkeys_min,skey_blocks,1);
    pkeys_max = repmat(pkeys_max,skey_blocks,1);

    pkeys_min = pkeys_min'; 
    pkeys_min = pkeys_min(:);    

    pkeys_max = pkeys_max';
    pkeys_max = pkeys_max(:);

    %skeys_min = skeys_min';
    skeys_min = skeys_min(:);

    %skeys_max = skeys_max';
    skeys_max = skeys_max(:);

    block_keys = [pkeys_min, pkeys_max, skeys_min, skeys_max];
    n_blocks = size(block_keys,1);
    
end
