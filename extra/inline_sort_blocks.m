function [] = inline_sort_blocks(job_meta_path)
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
%% FUnction : Makes a inline sorted block index and add to job meta file, currenly only the xl sorted one is available.

% Input: 
%       path of job meta file

% Writes to disk: 
%       Adds entry to job meta file : an array of live blocks which is inline sorted.


%%

job_meta = load(job_meta_path);                 % Load job meta information

block_keys = job_meta.block_keys;               % Block corner coordinates
block_keys =[ block_keys  (1:size(block_keys,1))' zeros(size(block_keys,1),1)]; % Append Block Keys S. No. at the end
[m n]=size(block_keys);

inc_xl=find(block_keys(:,3)>block_keys(1,3),1,'first')-1;                       % Find the xal-line increment
M=m/inc_xl;                                                                     % Find  the number of divisions along inline Currently this can also be found by squre rooting n
N=m/M;                                                                          % Find  the number of divisions along X-line
key_mat=zeros(M,N);                                                             
block_index = zeros(size(block_keys,1),2);                                      % Initialize a look-up table matrix

%%
% Loop for arranging the blocks in the look up table in a
for i=1:M
    for j=1:N
    key_mat(i,j)=j+N*(i-1);
    end
end
block_index(:,2)= key_mat(:);                                                   % Block index in x-line sorted format
key_mat_t=key_mat';                                                             % Transpose the lookup tabel
block_index(:,1)= key_mat_t(:);                                                 % Block index in in-line sorted format
job_meta.liveblocks_il_sorted=zeros(length(job_meta.liveblocks),1);             % Place for writing the blocks in an inline sorted fashion

%%
% Loop for sort thr live blocks in x-line sorted fashion
j=1;
k=-99;
for i=1:m
    k=find(job_meta.liveblocks==block_index(i,2));                              % find index of the live block(looping through the inline sorted blocks) in the x-line sorted live blocks rray
    %if you find the block
    if k > -99                                                                    
     job_meta.liveblocks_il_sorted(j) = job_meta.liveblocks(k) ;                % put it in the new array of in-line sorted live blocks
     j=j+1;                                                                     % increment counter
     k=-99;
    end  
end
block_keys2=job_meta.block_keys(job_meta.liveblocks_il_sorted,:);
save(job_meta_path,'-struct','job_meta','-v7.3');   % Saves Seismic structure to mat file
end
