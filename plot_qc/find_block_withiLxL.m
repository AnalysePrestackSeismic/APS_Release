function [ live_block_key ] = find_block_withiLxL( job_meta_path,iL,xL)
%UNTITLED Summary of this function goes here
%This find out the block that contains a user supplied Inline (iL) and xl
%(X-Line Number)
%if only using an inline use xL = 0
%if only using a xline use iL = 0
%   Detailed explanation goes here

iL = str2double (iL);
xL = str2double (xL);
job_meta = load(job_meta_path);



if(iL > 0 && xL>0)% While using both inline and xline
        flag_block_keys = job_meta.block_keys(job_meta.liveblocks,1)<iL & job_meta.block_keys(job_meta.liveblocks,2)>iL & job_meta.block_keys(job_meta.liveblocks,3)<xL & job_meta.block_keys(job_meta.liveblocks,4)>xL; % returns 1 if the il and xL is in the block (live block division scheme) else 0
        
elseif (xL == 0)% While using only inline
    flag_block_keys = job_meta.block_keys(job_meta.liveblocks,1)<iL & job_meta.block_keys(job_meta.liveblocks,2)>iL;
        
elseif (iL == 0)% While using only xline
    flag_block_keys = job_meta.block_keys(job_meta.liveblocks,3)<xL & job_meta.block_keys(job_meta.liveblocks,4)>xL;
    
end

live_block_key = find( flag_block_keys == 1); % Find the indexes associated with the blocks that has given Lines

live_block_key = num2str (live_block_key);% converting output to string

if size(live_block_key) == 0
    fprintf('No live block found with the given inline or xline'); % Telling user if the given inline or xline was not present in the live blocks
end

end

