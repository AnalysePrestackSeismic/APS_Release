function [ ] = segy_block_size( job_meta_path,varargin )
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
% segy_block_size: updates block sizes in job meta file
%
%   Arguments:
%       job_meta_path = pathname to job meta file
%
%   Optional arguments:
%       aperture = number of inlines and crosslines to overlap
%       inlines per block = force number of inlines per block
%       crosslines per block = force number of crosslines per block
%       (defaults to number of inlines per block)
%       


job_meta = load(job_meta_path);

%  check number of inlines/crosslines in volume

vol_ils = 1 + (job_meta.pkey_max - job_meta.pkey_min) / job_meta.pkey_inc;
vol_xls = 1 + (job_meta.skey_max - job_meta.skey_min) / job_meta.skey_inc;

% if block sizes are fixed then no need to call segy_make_blocks

if size(varargin,2)>=2
    block_ils = min(varargin{2},vol_ils);
    block_xls = min(varargin{end},vol_xls); % limit block size to vol size
    
    n_blocks_il = ceil(vol_ils/block_ils);
    n_blocks_xl = ceil(vol_xls/block_xls);
    job_meta.n_blocks = n_blocks_il * n_blocks_xl;
    
    il_blocks_min = job_meta.pkey_min + [0:n_blocks_il-1]'.*block_ils;
    il_blocks_max = il_blocks_min + block_ils - job_meta.pkey_inc;
    xl_blocks_min = job_meta.skey_min + [0:n_blocks_xl-1]'.*block_xls;
    xl_blocks_max = xl_blocks_min + block_xls - job_meta.skey_inc;
    
    job_meta.block_keys = zeros(job_meta.n_blocks,4);
    job_meta.block_keys(:,1) = sort(repmat(il_blocks_min,n_blocks_xl,1));
    job_meta.block_keys(:,2) = sort(repmat(il_blocks_max,n_blocks_xl,1));
    job_meta.block_keys(:,3) = repmat(xl_blocks_min,n_blocks_il,1);
    job_meta.block_keys(:,4) = repmat(xl_blocks_max,n_blocks_il,1);
   
else
    
    [job_meta.block_keys,job_meta.n_blocks] = segy_make_blocks(job_meta_path,varargin{:});
    
end


save(job_meta_path,'-struct','job_meta','-v7.3'); % Saves Seismic structure to mat file

% ##################################################################
% Find live blocks

% Assume that if blocks cover whole inline or crossline range then all are
% live

if block_ils == vol_ils || block_xls == vol_xls
    job_meta.liveblocks = [1:job_meta.n_blocks]';
else
    job_meta.liveblocks = select_live_blocks(job_meta_path);
end


save(job_meta_path,'-struct','job_meta','-v7.3'); % Saves Seismic structure to mat file

fprintf('Saved seismic structure to file ...\n')

end

