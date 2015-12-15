function segy_make_job(filepath,filename_string,il_byte,xl_byte,...
    offset_byte,parallel,output_dir,pkey_blocks,skey_blocks,tkey_blocks)
% -------------------------------------------------------------------------
% SEGY_MAKE_STRUCTURE: function to scan SEGY file to gain geometry and
% sample information used by other functions.
%   Inputs:
%       filepath = path of directory containing input angle stacks only as
%       cell array.
%       filename = name of SEGY file to scan.
%       il_byte  = inline number byte location
%       xl_byte  = crossline number byte location
%       output_dir = directory in which all DIGI outputs should be saved.
%   Outputs:
%       .mat file = metadata including sample rate, n_samples etc.
%       .mat_lite file = binary file containing IL/XL byte locations.
% -------------------------------------------------------------------------
pkey_loc = 1; % column numbers needs to be implemented
skey_loc = 2;
byte_loc = 3;
skey_max_loc = 4;
skey_inc_loc = 5;
tkey_loc = 6;
tkey_max_loc = 7;
tkey_inc_loc = 8;

il_byte = str2double(il_byte);
xl_byte = str2double(xl_byte);
offset_byte = str2double(offset_byte);
parallel = str2double(parallel);

[files_in,nfiles] = directory_scan(filepath,filename_string);
files_in.names = sort_nat(files_in.names);
start_point = pwd; % remember starting directory
filenames = '';

if parallel == 1
    %slurm submit
    % run using a little perl threading program
    for i_file = 1:1:nfiles
        filenames = horzcat(filenames,' ',files_in.names{i_file});
    end
    filepath = files_in.path{1};
    [result, status] = perl('/apps/gsc/scripts/perl_matlab_par_for.pl','/apps/gsc/matlab-mcode-beta/eslib/final_digi_condensed/scan_segy/run_segy_make_structure.sh','/apps/matlab/v2011a',filepath,'189','193',filenames);
    
else
    for i_file = 1:1:nfiles
        
        filename = files_in.names{i_file};
        filepath = files_in.path{i_file};
        segy_make_structure(filepath,num2str(il_byte),num2str(xl_byte),filename);
    end
    
end


%job_meta.non_live_traces = non_live_traces;
% Save meta information about files scanned to .mat file
for i_file = 1:1:size(files_in.names,2)
    if ~isempty(strfind(files_in.names{i_file},'segy'))
        job_meta.files{i_file} = regexprep(files_in.names{i_file}, 'segy', 'mat_orig_lite');
    elseif ~isempty(strfind(files_in.names{i_file},'sgy'))
        job_meta.files{i_file} = regexprep(files_in.names{i_file}, 'sgy', 'mat_orig_lite');
    end
end
%job_meta.files = cell2mat(job_meta.files);
job_meta.files = unique(job_meta.files);
job_meta.paths = unique(files_in.path');
job_meta.output_dir = output_dir;
job_meta.il_byte = il_byte;
job_meta.xl_byte = xl_byte;
job_meta.offset_byte = offset_byte;

vol_names = strfind(files_in.names', '_block');

for i_file = 1:1:nfiles
    job_meta.volumes{i_file} = files_in.names{i_file}(1:vol_names{i_file}-1);
end
job_meta.volumes = unique(job_meta.volumes)';
job_meta.nvols = size(job_meta.volumes,1);
for i_vol = 1:1:job_meta.nvols
    job_meta.angle{i_vol,1} = str2double(regexp(job_meta.volumes{i_vol},'(\d{2})','match'));
    
    % find files associated with volume from all blocks
    [files_in,nfiles] = directory_scan(job_meta.paths,job_meta.volumes{i_vol});
    job_meta.vol_traces{i_vol,1} = 0;
    ii = 1;
    % loop round all the mat_lite files to do with this volume
    for il = 1:nfiles
        if strfind(files_in.names{il},'.mat_orig_lite')
            seismic = segy_read_binary(strcat(files_in.path{il},files_in.names{il}));
            if ii == 1
                vol_index{i_vol} = seismic.trace_ilxl_bytes;
            else
                vol_index{i_vol} = [vol_index{i_vol}; seismic.trace_ilxl_bytes];
            end
            
            job_meta.vol_traces{i_vol,1} = job_meta.vol_traces{i_vol,1}+seismic.n_traces;
            pkey_min(ii) = min(seismic.trace_ilxl_bytes(:,pkey_loc));
            pkey_max(ii) = max(seismic.trace_ilxl_bytes(:,pkey_loc));
            pkey_inc(ii) = mode(diff(unique(seismic.trace_ilxl_bytes(:,pkey_loc))));
            skey_min(ii) = min(seismic.trace_ilxl_bytes(:,skey_loc));
            skey_max(ii) = max(seismic.trace_ilxl_bytes(:,skey_max_loc));
            skey_inc(ii) = mode(seismic.trace_ilxl_bytes(:,skey_inc_loc));
            
            if seismic.is_gather == 1
                file_pkey_start{i_vol}(ii,:) = [ seismic.trace_ilxl_bytes(1,pkey_loc) seismic.trace_ilxl_bytes(end,pkey_loc) seismic.trace_ilxl_bytes(1,skey_loc) seismic.trace_ilxl_bytes(end,skey_max_loc) seismic.trace_ilxl_bytes(1,tkey_loc) seismic.trace_ilxl_bytes(end,tkey_loc) uint64(files_in.names{il})];
            else
                file_pkey_start{i_vol}(ii,:) = [ seismic.trace_ilxl_bytes(1,pkey_loc) seismic.trace_ilxl_bytes(end,pkey_loc) seismic.trace_ilxl_bytes(1,skey_loc) seismic.trace_ilxl_bytes(end,skey_max_loc) uint64(files_in.names{il}) ];
            end
            % read filename with char(file_pkey_start{1,ivol}(ii,5:end))
            
            if seismic.is_gather == 1
                tkey_min(ii) = min(seismic.trace_ilxl_bytes(:,tkey_loc));
                tkey_max(ii) = max(seismic.trace_ilxl_bytes(:,tkey_max_loc));
                tkey_inc(ii) = mode(seismic.trace_ilxl_bytes(:,tkey_inc_loc));
            end
            ii = ii + 1;
        end
    end
    job_meta.n_samples{i_vol} = seismic.n_samples;
    job_meta.trc_head{i_vol} = 240;
    job_meta.bytes_per_sample{i_vol} = 4;
    job_meta.vol_nblocks(i_vol,1) = ii-1;
    job_meta.pkey_min(i_vol,1) = min(pkey_min);
    job_meta.pkey_max(i_vol,1) = max(pkey_max);
    job_meta.pkey_inc(i_vol,1) = mode(pkey_inc);
    job_meta.skey_min(i_vol,1) = min(skey_min);
    job_meta.skey_max(i_vol,1) = max(skey_max);
    job_meta.skey_inc(i_vol,1) = mode(skey_inc);
    if seismic.is_gather == 1
        job_meta.tkey_min(i_vol,1) = min(tkey_min);
        job_meta.tkey_max(i_vol,1) = max(tkey_max);
        job_meta.tkey_inc(i_vol,1) = mode(tkey_inc);
        job_meta.is_gather = 1;
    else
        job_meta.is_gather = 0;
    end
end

% ############# code to find missing traces and remake mat lite files

% now loop through each pkey from pkey_min in pkey_inc to find only the
% matching skey locations

%pkeys = unique(seismic.trace_ilxl_bytes(:,1));
pkeyn = 1+((mode(job_meta.pkey_max)-mode(job_meta.pkey_min))/mode(job_meta.pkey_inc));
skeyn = 1+((mode(job_meta.skey_max)-mode(job_meta.skey_min))/mode(job_meta.skey_inc));
pkey_min = mode(job_meta.pkey_min);
skey_min = mode(job_meta.skey_min);
pkey_max = mode(job_meta.pkey_max);
skey_max = mode(job_meta.skey_max);
pkey_inc_mode = mode(job_meta.pkey_inc);
skey_inc_mode = mode(job_meta.skey_inc);
pkey_log_idx = false(pkeyn,job_meta.nvols);
% find all the unique pkeys
for i_vol = 1:1:job_meta.nvols
    for i_pkey = 1:1:size(vol_index{i_vol},1)
        pkey_log_idx(((vol_index{i_vol}(i_pkey,pkey_loc) - pkey_min)/pkey_inc_mode + 1),i_vol) = 1;
    end
end
pkey_vals(:,1) = pkey_min:pkey_inc_mode:pkey_max;
pkey_vals_unif = pkey_vals(sum(pkey_log_idx,2) == job_meta.nvols,:);

% loop round for each pkey to find the unique skeys
skey_vals(:,1) = skey_min:mode(job_meta.skey_inc):skey_max;
count = 1;
for i_pkey = 1:1:size(pkey_vals_unif,1)
    skey_log_idx = false(skeyn,job_meta.nvols);
    for i_vol = 1:1:job_meta.nvols
        il_rows = find(vol_index{i_vol}(:,1) == pkey_vals_unif(i_pkey,1));
        for s_entry = 1:1:size(il_rows,1)
            skey_start = (vol_index{i_vol}(il_rows(s_entry),skey_loc)-skey_min)/skey_inc_mode+1;
            skey_end = (vol_index{i_vol}(il_rows(s_entry),skey_max_loc)-skey_min)/skey_inc_mode+1;
            %skey_inc = vol_index{i_vol}(il_rows(s_entry),skey_inc_loc);
            skey_log_idx(skey_start:skey_end,i_vol) = 1;
            % this was edited by JS on 02/12/2013 based on Uruguay data
            % (SW)
        end
    end
    skey_vals_unif = skey_vals(sum(skey_log_idx,2) == job_meta.nvols,:);
    
    % Compress skey_vals_unif on the fly
    start_idx = 1;
    blocktr = size(skey_vals_unif,1);
    row_i = 0;
    while start_idx < blocktr
        cdp_1 = skey_vals_unif(start_idx);
        if start_idx < blocktr-1
            cdp_2 = skey_vals_unif(start_idx+1);
        end
        if (cdp_2 - cdp_1) == skey_inc_mode
            row_i = row_i+1;
        else
            compress_ilxl_bytes(count,1) = pkey_vals_unif(i_pkey,1);
            compress_ilxl_bytes(count,2) = skey_vals_unif(start_idx-row_i,1);
            compress_ilxl_bytes(count,3) = 0; % byte locations later
            
            if cdp_2 == cdp_1
                compress_ilxl_bytes(count,4) = skey_vals_unif(start_idx+1,1);
            else
                compress_ilxl_bytes(count,4) = skey_vals_unif(start_idx,1);
            end
            
            xl_inc = skey_vals_unif(start_idx-row_i+1,1)-skey_vals_unif(start_idx-row_i,1);
            compress_ilxl_bytes(count,5) = xl_inc;
            
            row_i = 0;
            count = count + 1;
        end
        start_idx = start_idx + 1;
    end
end

% loop over each volume and output a compressed mat_lite file for each sub
% file
for i_vol = 1:1:job_meta.nvols
    trc_length =  (job_meta.n_samples{i_vol}*job_meta.bytes_per_sample{i_vol})+...
        job_meta.trc_head{i_vol};
    
    % loop over all the sub files in each volume
    %mlf_pkey_start = file_pkey_start{1,ivol}(file_loop,1);
    %file_name = char(file_pkey_start{1,ivol}(file_loop,5:end));
    file_loop = 1;
    
    for i_key = 1:1:size(compress_ilxl_bytes,1)
        
        % loop round all entries in the files table for this volume and
        % find which file the pkey is in and then store the file name and
        % compare when it changes
        
        for file_loop = 1:1:size(file_pkey_start{1,i_vol},1)
            if compress_ilxl_bytes(i_key,pkey_loc) > file_pkey_start{1,i_vol}(file_loop,1) && compress_ilxl_bytes(i_key,pkey_loc) < file_pkey_start{1,i_vol}(file_loop,2);
                cur_output_filename{i_key} = [char(file_pkey_start{1,i_vol}(file_loop,5:end))];
            elseif compress_ilxl_bytes(i_key,pkey_loc) == file_pkey_start{1,i_vol}(file_loop,1) || compress_ilxl_bytes(i_key,pkey_loc) == file_pkey_start{1,i_vol}(file_loop,2);
                
                % need to test for the skey being in the range of the file
                % or maybe the next file
                test_cur_skey_min = compress_ilxl_bytes(i_key,skey_loc);
                test_cur_skey_max = compress_ilxl_bytes(i_key,skey_max_loc);
                if file_pkey_start{1,i_vol}(file_loop,3) >= test_cur_skey_min && file_pkey_start{1,i_vol}(file_loop,4) <= test_cur_skey_max
                    cur_output_filename{i_key} = [char(file_pkey_start{1,i_vol}(file_loop,5:end))];
                else
                    cur_output_filename{i_key} = 'hfysorth';
                end
            else
                cur_output_filename{i_key} = 'hfysorth';
            end
        end
        
        il_rows = find(vol_index{i_vol}(:,pkey_loc) == compress_ilxl_bytes(i_key,pkey_loc));
        
        il_found = vol_index{i_vol}(il_rows,:);
        
        xl_rows = (il_found(:,2) <= compress_ilxl_bytes(i_key,2)) & (il_found(:,4) >= compress_ilxl_bytes(i_key,2));
        
        xl_found = il_found(xl_rows,:);
        
        [xl,ind] = min(xl_found(:,2));
        check_xl = (compress_ilxl_bytes(i_key,2) - xl)/xl_found(ind,5);
        
        n_traces_away = check_xl;
        n_bytes_away = n_traces_away*(trc_length);
        
        compress_ilxl_bytes(i_key,3) = xl_found(ind,3)+n_bytes_away;
        
    end
    
    
    % write out the new compressed files
    
    prev_fileout = 'jsdjkfhjsiso';
    i_filep = 0;
    
    for i_key = 1:1:size(compress_ilxl_bytes,1)
        % write out the new compressed files
        if strcmp('hfysorth',cur_output_filename{i_key}) == 0  % strings are not equal, case sensitive)
            if strcmp(prev_fileout,cur_output_filename{i_key}) == 0  % strings are not equal, case sensitive)
                i_filep = i_filep + 1;
                byte_loc2{i_filep,1}(i_key,:) = compress_ilxl_bytes(i_key,:);
                byte_loc2{i_filep,2} = cur_output_filename{i_key};
                prev_fileout = cur_output_filename{i_key};
            else
                byte_loc2{i_filep,1}(i_key,:) = compress_ilxl_bytes(i_key,:);
            end
        end
    end
    
    for f_key = 1:1:size(byte_loc2,1)
        %loop and write out all the file mat_lite files
        
        % need to index to the seismic file name correctly  ##########
        % and cannot get the fopen to work, works on command line though
        seismic = segy_read_binary(strcat(job_meta.paths{1},job_meta.files{i_vol}));
        filepath_binary = uint64(seismic.filepath);
        pad_filepath = zeros(1,(2000-length(filepath_binary)));
        filepath_binary = [filepath_binary,pad_filepath];
        write_file = strcat(job_meta.paths{1},strrep(job_meta.files{i_vol},'.mat_orig_lite','.mat_lite'));
        fid_writecj = fopen(write_file,'w');
        
        fwrite(fid_writecj,[filepath_binary';seismic.file_type;seismic.s_rate;seismic.n_samples;seismic.n_traces;seismic.pkey;seismic.skey;seismic.tkey;seismic.is_gather],'double');
        
        fwrite(fid_writecj,reshape(byte_loc2{f_key,1}(:,:)',[],1),'double');
        fclose(fid_writecj);
    end
    
end
% do we need to do any of this for single volumes!
% ##################################################################
job_meta.files = regexprep(job_meta.files,'mat_orig_lite','mat_lite')';
%job_meta.files = reshape(job_meta.files,[],job_meta.nvols);
job_meta.s_rate = seismic.s_rate;
job_meta.vol_traces = cell2mat(job_meta.vol_traces);
% store mean angle?
% separate blocks for a single volume
%job_meta.volumes = cell2mat(job_meta.volumes);
str_date = date;
str_date = regexprep(str_date, '-', '');
job_meta_dir = strcat(job_meta.output_dir,'job_meta/');
mkdir(job_meta_dir);
job_meta_path = strcat(job_meta_dir,'job_meta_',str_date,'.mat');

save(job_meta_path,'-struct','job_meta','-v7.3'); % Saves Seismic structure to mat file
job_meta = load(job_meta_path);
[job_meta.block_keys,job_meta.n_blocks] = segy_make_blocks(job_meta_path,pkey_blocks,skey_blocks,tkey_blocks);

save(job_meta_path,'-struct','job_meta','-v7.3'); % Saves Seismic structure to mat file

fprintf('Saved seismic structure to file ...\n')

cd(start_point)
end