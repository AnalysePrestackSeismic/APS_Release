function [] = select_live_blocks(job_meta_path)
%SELECT_LIVE_BLOCKS Summary of this function goes here
%   this reads all the blocks and makes an array of live iblocks where live
%   is a number of traces live greater than 10, this takes the middle
%   volume to look for data

job_meta = load(job_meta_path);
%iloop = 1;
%liveb = zeros(job_meta.n_blocks,1);
% Read traces for 1/2 max angle stack
%vol_index_wb = ceil(job_meta.nvols*0.5);

% need to loop on iblock
col_count = 5;
max_file_number = -100;
if job_meta.is_gather == 1
    for i_block = 1:1:size(job_meta.files,2)
        seismic_mat_path = [job_meta.paths{1}, job_meta.files{i_block}];
        seismic = segy_read_binary(seismic_mat_path);
        
        min_il = min(seismic.trace_ilxl_bytes(:,1));
        max_il = max(seismic.trace_ilxl_bytes(:,1));
        min_xl = min(seismic.trace_ilxl_bytes(:,2));
        max_xl = max(seismic.trace_ilxl_bytes(:,4));
        
        file_number_str = strfind(job_meta.files{i_block},'k');
        file_number_end = strfind(job_meta.files{i_block},'.');
        file_number = str2num(job_meta.files{i_block}(file_number_str+1:file_number_end-1));
        
        rows_log = job_meta.block_keys(:,1) >= min_il ...
            & job_meta.block_keys(:,2) <= max_il & ...
            job_meta.block_keys(:,3) >= min_xl & ...
            job_meta.block_keys(:,4) <= max_xl;
        
        job_meta.block_keys(rows_log,col_count) = file_number;
        col_count = col_count + 1;
        
        if max_file_number < file_number
            max_file_number = file_number;
        end
    end
    
    % find locations that have one file entry
    block_files = sum(job_meta.block_keys(sum(job_meta.block_keys(:,5:end),2) ...
        <= max_file_number,5:end),2);
    
    % add cols for locations with multi-files
    
    
else
    
    i_counter = 1;
    
    for i_block = 1:1:size(job_meta.files,1)
        seismic_mat_path = [job_meta.paths{1}, job_meta.files{i_block}];
        seismic = segy_read_binary(seismic_mat_path);
        expand_keys = job_meta.block_keys;
        %
        for i_key = 1:1:size(expand_keys,1)
            vol_keys = segy_index_byte_finder(job_meta_path,expand_keys(i_key,:),'1');
            
            if size
            job_meta.liveblocks(i_counter) = i_key;
            i_counter = i_counter + 1;
%             il_rows = seismic.trace_ilxl_bytes(:,1) == expand_keys(i_key,1);
%             
%             if sum(il_rows) > 0
%                 il_found = seismic.trace_ilxl_bytes(il_rows,:);
%                 xl_rows = (il_found(:,2) <= expand_keys(i_key,3)) & (il_found(:,4) >= expand_keys(i_key,4));
%                 
%                 if sum(xl_rows) > 0
%                     % found an skey     
%                     job_meta.liveblocks(i_counter) = i_key;
%                     i_counter = i_counter + 1;
%                 else
%                     % no skey found
%                     
%                 end
%             else
%                 % no peky found
%             end
         end
     end
    
end
        
%         seismic_mat_path = [job_meta.paths{1}, job_meta.files{i_block}];
%         seismic = segy_read_binary(seismic_mat_path);
%         
%         for i_block_row = 1:1:size(job_meta.block_keys,1)
%             min_il = job_meta.block_keys(i_block_row,1);
%             max_il = job_meta.block_keys(i_block_row,2);
%             min_xl = job_meta.block_keys(i_block_row,3);
%             max_xl = job_meta.block_keys(i_block_row,4);
%             
%             rows_log
%             
%         end
%         
%         min_il = min(seismic.trace_ilxl_bytes(:,1));
%         max_il = max(seismic.trace_ilxl_bytes(:,1));
%         min_xl = min(seismic.trace_ilxl_bytes(:,2));
%         max_xl = max(seismic.trace_ilxl_bytes(:,4));
%         
%         rows_log = job_meta.block_keys(:,1) >= min_il ...
%             & job_meta.block_keys(:,2) <= max_il & ...
%             job_meta.block_keys(:,3) >= min_xl & ...
%             job_meta.block_keys(:,4) <= max_xl;
%         
%         block_numbers = (1:1:size(job_meta.block_keys,1))';
%         live_blocks(:,i_block) = block_numbers(rows_log);
        
        %[~, traces,~,~] = node_segy_read(job_meta_path,num2str(vol_index_wb),num2str(i_block));
        
        %[n_samples,n_traces] = size(traces);
        
        %if n_traces > 10
        % write out an accepetable block number into the output array
        %    liveb(iloop) = i_block;
        %iloop = iloop + 1;
        %end

    parallel = 0;
    if parallel == 1
        % run using a little perl threading program
        % need to compile the wrpper round the node_segy_read code that
        % does the size and return the iblock if it is non blank otherwise
        % return zero.
        % then sort all the replies and drop the zeros
        for i_file = 1:1:nfiles
            filenames = horzcat(filenames,' ',files_in.names{i_file});
        end
        filepath = files_in.path{1};
        [result, status] = perl('/apps/gsc/scripts/perl_matlab_par_for.pl','/apps/gsc/matlab-library/final_digi_condensed/scan_segy/run_segy_make_structure.sh','/apps/matlab/v2011a',filepath,'189','193',filenames);
    end

% truncate the live block array before writting

%job_meta.liveblocks = liveb(1:(iloop-1));

save(job_meta_path,'-struct','job_meta','-v7.3'); % Saves Seismic structure to mat file
%find_job_name = strfind(job_meta_path,'/');
% write out the array to file to use again
%save(strcat(job_meta.output_dir,[job_meta_path(find_job_name(end)+1:end-4) '_live_iblocks.mat']),'liveblocks','-v7.3');
%save(strcat(job_meta.output_dir,'live_iblocks.bin'),liveb(1:(iloop-1)));
%
% fid_wav = fopen(strcat(job_meta.output_dir,[job_meta_path(find_job_name(end)+1:end-4) '_live_iblocks.bin']),'w');
% fwrite(fid_wav,liveb(1:(iloop-1)),'uint32');
% fclose(fid_wav);

end

