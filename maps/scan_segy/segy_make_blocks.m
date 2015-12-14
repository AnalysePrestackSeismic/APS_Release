function [block_keys,n_blocks] = segy_make_blocks(job_meta_path,varargin)
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
%   segy_make_blocks: makes _mat_lite and mat_orig_lite files from the Seismic Header 
%   Run this on one angle stack as they should all have the same
%   geometry/file structure. This will create a mat file in the location of the input SEGY definining the
%   structure of the input file.Structure format:
%   PKey SKey Byte_Loc SKey_max SKey_inc TKey TKey_max TKey_inc
%   Use segy_index_checker(seismic_mat_path) to check scan validity
%   
%   Arguments:  
%       filepath = path of segy file
%       il_byte = inline byte location
%       xl_byte = Cross Line byte Location
%       offset_byte =  Offset Byte Location
%       anggath = 1 for angle gathers and anything else for offset gathers
%       shots might have a different number
%       filename = Name of Segy File
%   Optional arguments:
%       aperture = number of inlines and crosslines to overlap
%       inlines per block = force number of inlines per block
%       crosslines per block = force number of crosslines per block
%       (defaults to number of inlines per block)
%
%   Note: Typical IL/XL byte locations: 189 & 193.
%   Offset byte location is typically 37 
%
%   Outputs:
%       .mat file = metadata including sample rate, n_samples etc.
%       .mat_lite file = binary file containing IL/XL byte locations.
%
%   Writes to Disk:
%       job meta files: give description and paths

%%
    job_meta = load(job_meta_path);
    
    
    aperture=0; pkey_inc=0; skey_inc=0;
    
    if size(varargin,2) >= 1 
        aperture = varargin{1};
    end
    
    if size(varargin,2) >= 2 
        pkey_inc = varargin{2};
        skey_inc = pkey_inc;
    end
    
    if size(varargin,2) >= 3
        skey_inc = varargin{3};
    end
    
    if pkey_inc == 0
        % this is the size of memory that is defined per block, needs to make
        % sure that the n cores per node * gigabytes_per_block is less than the
        % memory per node
        
        gigabytes_per_block = 2;
        
        total_blocks = ceil(sum(job_meta.vol_traces.*((cell2mat(job_meta.n_samples)'*4)+240)/1024/1024/1024)/gigabytes_per_block);
        
        if total_blocks < 1104;
            if total_blocks < 1000;
                total_blocks = 1006;
            else
                total_blocks = 1102;
            end
        end
        
        
        %work out if this means that there are enough traces in a block
        
        pkey_blocks = ceil(sqrt(total_blocks));
        skey_blocks = pkey_blocks;
        %tkey_blocks = 1;
        
        % Need to correct
        %pkey_inc = floor((mode(job_meta.pkey_max) - mode(job_meta.pkey_min))/pkey_blocks);
        pkey_inc = floor(((mode(job_meta.pkey_max) - mode(job_meta.pkey_min))/pkey_blocks)/job_meta.pkey_inc)*job_meta.pkey_inc;
        %skey_inc = floor((mode(job_meta.skey_max) - mode(job_meta.skey_min))/skey_blocks);
        skey_inc = floor(((mode(job_meta.skey_max) - mode(job_meta.skey_min))/skey_blocks)/job_meta.skey_inc)*job_meta.skey_inc;
        
    end
    
    if pkey_inc == 0
        pkeys_min = min(job_meta.pkey_min);
    else
        pkeys_min = min(job_meta.pkey_min):pkey_inc:max(job_meta.pkey_max)-pkey_inc;
    end
    
    if skey_inc == 0
        skeys_min = min(job_meta.skey_min);
    else
        skeys_min = min(job_meta.skey_min):skey_inc:max(job_meta.skey_max)-skey_inc;
    end
    
    pkey_blocks = size(pkeys_min,2);
    skey_blocks = size(skeys_min,2);
    
    pkeys_max = pkeys_min - mode(job_meta.pkey_inc);
    pkeys_max = pkeys_max(2:end);
    pkeys_max(end+1) = mode(job_meta.pkey_max);

    skeys_max = skeys_min - mode(job_meta.skey_inc);
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

    block_keys = [pkeys_min-aperture, pkeys_max+aperture, skeys_min-aperture, skeys_max+aperture];
    n_blocks = size(block_keys,1);    
    
end
