function job_meta_path = segy_make_job_rg(filepath,filename_string,il_byte,xl_byte,...
    offset_byte,parallel,anggath,output_dir)
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
% segy_make_job: function to scan SEGY file to gain geometry and
% sample information used by other functions.

% File Name Convention:
% Angle Stacks: startangle-endangle_name_block1.sgy where 1 is the volume number and would be incremented if there are multiple segy files per volume, for angles <10, use format 01,02,03 etc,e.g. 02_05 not 2_5.
% Gather: gather_name_block1.sgy where 1 is the volume number and would be incremented if there are multiple segy files per volume
% File names shouldn't have substring 'segy' since its an identifier in the program 
%
%   Arguments:
%       filepath =          path of directory containing input angle stacks only as cell array, {'/path/'}
%       filename_string =   name of search-string in SEGY file name to scan
%       il_byte  =          inline number byte location (189 for SEGY Rev1 format)
%       xl_byte  =          crossline number byte location (193 for SEGY Rev1 format)
%       offset_byte =       Offset byte location (37 for SEGY Rev1 format)
%       parallel =          1 if you want to run in parallel O if you want to run in a single machine
%       anggath =           1 if angle gathers
%       output_dir =        directory in which all DIGI outputs should be saved.
%   
%   Outputs:
%       .mat file = metadata including sample rate, n_samples etc.
%       .mat_lite file = binary file containing IL/XL byte locations.
%
%   Writes to Disk:
%       job meta files: give description and paths

% Note: Search for the string 'gather' in the file name.

%%
%-------------------PROCESSING FUNCTION ARGUMENTS--------------------------
% Column numbers define output format of .mat_lite file. 
% Should make this a global format definition.  
pkey_loc = 1;                               % column numbers needs to be implemented Primary Key
skey_loc = 2;                               % Secondary Key
byte_loc = 3;                               % Byte location
skey_max_loc = 4;                           % Secondary Key Maximum
skey_inc_loc = 5;                           % Secondary Key Increment  
tkey_loc = 6;                               % Tertiary Key
tkey_max_loc = 7;                           % Tertiary Key Maximum
tkey_inc_loc = 8;                           % Tertiary Key Increment 

il_byte = str2double(il_byte);              % catching and converting to number inline byte location from function argument
xl_byte = str2double(xl_byte);              % catching and converting to number xline byte location from function argument     
offset_byte = str2double(offset_byte);      % catching and converting to number offset byte location
parallel = str2double(parallel);            % catching and converting to number to find whether to run function in parallel or not
anggath = str2double(anggath);

if ~iscell(filepath)
    filepath_cell{1} = filepath;
    clear filepath
    filepath = filepath_cell;
    clear filepath_cell;
end
[files_in,nfiles] = directory_scan(filepath,filename_string); % files_in is a structure, .names cell array, .path cell array 
                                                              % directory_scan filters out the file names with the user supplied string and also returns the number of such files: nfiles that ...  
                                                              % ... can be used for index pre alocation in future steps 
files_in.names = sort_nat(files_in.names);                    % Natural Sort file names 
start_point = pwd;                                            % Remember starting directory
%%                
%-------------------SCAN THROUGH FILTERED SEGY FILES ----------------------
%------------------CREATE .MAT_ORIG_LITE FILE AFTER SCANNING SEGY--------
filenames = '';                                                     % Initialize Filenames
if parallel == 1                                                    % For running in parallel
    %slurm submit: run using a perl threading program
    for i_file = 1:1:nfiles
        filenames = horzcat(filenames,' ',files_in.names{i_file});  % Create a single string of all filtered filenames
    end
    filepath = files_in.path{1};
    [result, status] = perl('/apps/gsc/matlab-library/development/maps/extra/perl_matlab_par_for.pl','/apps/gsc/matlab-library/development/maps/scan_segy/run_segy_make_structure.sh','/apps/matlab/v2013a',filepath,num2str(il_byte),num2str(xl_byte),num2str(offset_byte),num2str(anggath),filenames);
    
else                                                                % For running in single machine
    for i_file = 1:1:nfiles
        
        filename = files_in.names{i_file};                          % Sequentially put file name in the structure
        filepath = files_in.path{i_file};                           % Sequentially put file path in the structure
        segy_make_structure(filepath,num2str(il_byte),num2str(xl_byte),num2str(offset_byte),num2str(anggath),filename); % Scan segy and make mat file <file name.mat_orig_lite> with structure ...
                                                                                                       %... of format [ PKey SKey Byte_Loc SKey_max SKey_inc TKey TKey_max TKey_inc ]...
                                                                                                       %... This checks if the file already exists. So if program is run repeatedly...
                                                                                                       %...it checks lif the scan already exists skips rescanning. This also makes a <file name.mat.lite >...
                                                                                                       %...file for later use
    end
    
end
%------------FIND WHETHER FILE IS A GATHER OR NOT----------------------
if ~isempty(strfind(files_in.names{i_file},'gath'))                 % Search for the string 'gath' in the file name
    is_gather = 1;                                                  % Check if file name says its a gather
else
    is_gather = 0;                                                  % Check if file name says its not a gather    
end
%-----------------------WRITE FILE NAMES IN JOB_META FILE---------------
%job_meta.non_live_traces = non_live_traces;
% Save meta information about files scanned to .mat file
for i_file = 1:1:size(files_in.names,2)                                                             % Scan through the filtered files sequentially. Number of files known from size of structure
    if ~isempty(strfind(files_in.names{i_file},'segy'))                                             % Check whether file name has segy in it
%        if is_gather == 0
            job_meta.files{i_file} = regexprep(files_in.names{i_file}, 'segy', 'mat_orig_lite');    % Replace 'segy' by 'mat_orig_lite' and puts in job_meta file
            if_string = 'mat_orig_lite';                                                            % String and Flag for later use
%         else
%             job_meta.files{i_file} = regexprep(files_in.names{i_file}, 'segy', 'mat_lite');
%             if_string = 'mat_lite';
%         end
    elseif ~isempty(strfind(files_in.names{i_file},'sgy'))                                          % Check whether file name has segy in it
%        if is_gather == 0
            job_meta.files{i_file} = regexprep(files_in.names{i_file}, 'sgy', 'mat_orig_lite');     % Replace 'segy' by 'mat_orig_lite' and put in job_meta file
            if_string = 'mat_orig_lite';                                                            % String and Flag for later use
%         else
%             job_meta.files{i_file} = regexprep(files_in.names{i_file}, 'sgy', 'mat_lite');
%             if_string = 'mat_lite';
%         end
    end
    
%     if ~isempty(strfind(files_in.names{i_file},'gath'))
%         is_gather = 1;
%     else
%         is_gather = 0;        
%     end
end
%%
%--------RESTRUCTURE JOB META FILE (REMOVE NON-ENTRIES, DUPLICATIONS, ETC) AND ENTER MORE INFO IN JOB_META FILE--------
%job_meta.files = cell2mat(job_meta.files);
count = 1;                                          % Initializing counter
for i_file = 1:1:size(job_meta.files,2)             % Loop for counting actual number of files filterd into job_meta file 
    if ~isempty(job_meta.files{i_file})
        files_tmp{count} = job_meta.files{i_file};  % Put non null entries in file_tmp
        count = count + 1;                          % Increment counter
    end     
end
job_meta.files = unique(files_tmp);                 % Remove duplicate entries (file names)  if any in job_meta file
job_meta.paths = unique(files_in.path');            % Remove duplicate entries in path names if any in job meta file
job_meta.output_dir = output_dir;                   % Write output directory from argument in job meta file
job_meta.il_byte = il_byte;                         % Write inline byte location from argument in job meta file
job_meta.xl_byte = xl_byte;                         % Write X-line byte from argument in job meta file
job_meta.offset_byte = offset_byte;                 % Write offset byte from argument in job meta file

vol_names = strfind(files_in.names', '_block');     %Find string '_block' in filtered file names and returns the index at which it starts in the file name

for i_file = 1:1:nfiles
    job_meta.volumes{i_file} = files_in.names{i_file}(1:vol_names{i_file}-1); % Enter the file name till the index vol_names into job_meta file  
end
job_meta.volumes = unique(job_meta.volumes)';       % Removes duplicate entries. Also helps if the function is run multiple times to filter out a .mat_orig_lite file with same file name
job_meta.nvols = size(job_meta.volumes,1);          % Finds and enters the total number of Volumes
%%
%----READ THE SEISMIC HEADER INFORMATION FROM  FILES----------

for i_vol = 1:1:job_meta.nvols                      % Loop run for all Volumes sequentially
    if is_gather == 1                               % If this is a gather file do nothing since there shouldnt be multiple gathers in one block
    else                                            % If this is an angle stack
        job_meta.angle{i_vol,1} = str2double(regexp(job_meta.volumes{i_vol},'(\d{2})','match'));    % Find the angle of the angle stack volume ? and put it in job_meta file
    end
    
    [files_in,nfiles] = directory_scan(job_meta.paths,job_meta.volumes{i_vol});                     % Find files associated with volume from all blocks            
    job_meta.vol_traces{i_vol,1} = 0;                                                               % Intialize number of trace ? to zero
    ii = 1;                                                                                         % Intialize index for following loop .
    
    %---loop round all the mat_lite files to do with this volume-------
    for il = 1:nfiles                                                                               % nfiles is no. of files retuned from scanning directory
        if strfind(files_in.names{il},if_string)                                                    % If a segy file was found this must have been initialized to '' else will be null?
            seismic = segy_read_binary(strcat(files_in.path{il},files_in.names{il}));               % Read seismic data related information from <file name. mat_orig_liet file> into a structure
            if ii == 1
                vol_index{i_vol} = seismic.trace_ilxl_bytes;
            else
                vol_index{i_vol} = [vol_index{i_vol}; seismic.trace_ilxl_bytes];
            end
            
            job_meta.vol_traces{i_vol,1} = job_meta.vol_traces{i_vol,1}+seismic.n_traces;
            
            % temp fix to fix a ury dataset
%             if min(seismic.trace_ilxl_bytes(:,skey_loc)) == 0
%                 min(seismic.trace_ilxl_bytes(:,skey_loc))
%                 strcat(files_in.path{il},files_in.names{il})
%             end
           
            %pkey_min(ii) = min(seismic.trace_ilxl_bytes(:,pkey_loc));
            pkey_min(ii) = min(seismic.trace_ilxl_bytes(seismic.trace_ilxl_bytes(:,pkey_loc) > 0,pkey_loc));
            pkey_max(ii) = max(seismic.trace_ilxl_bytes(:,pkey_loc));
            if pkey_min(ii) == pkey_max(ii)
                pkey_inc(ii) = 1;
            else
                pkey_inc(ii) = mode(diff(unique(seismic.trace_ilxl_bytes(:,pkey_loc))));
            end
            %skey_min(ii) = min(seismic.trace_ilxl_bytes(:,skey_loc));
            skey_min(ii) = min(seismic.trace_ilxl_bytes(seismic.trace_ilxl_bytes(:,skey_loc) > 0,skey_loc));
            skey_max(ii) = max(seismic.trace_ilxl_bytes(:,skey_max_loc));
            

            if skey_min(ii) == skey_max(ii)
                skey_inc(ii) = 1;
            else
                %skey_inc(ii) = mode(seismic.trace_ilxl_bytes(:,skey_inc_loc));
                skey_inc(ii) = mode(seismic.trace_ilxl_bytes((seismic.trace_ilxl_bytes(:,skey_max_loc) - seismic.trace_ilxl_bytes(:,skey_loc)) .* seismic.trace_ilxl_bytes(:,skey_inc_loc) > 0,skey_inc_loc));
            end
            file_pkey_start{i_vol}(ii,:) = zeros(1,100);
            
            if seismic.is_gather == 1
                length_il = size([ seismic.trace_ilxl_bytes(1,pkey_loc) seismic.trace_ilxl_bytes(end,pkey_loc) seismic.trace_ilxl_bytes(1,skey_loc) seismic.trace_ilxl_bytes(end,skey_max_loc) seismic.trace_ilxl_bytes(1,tkey_loc) seismic.trace_ilxl_bytes(end,tkey_loc) uint64(files_in.names{il})],2);
                file_pkey_start{i_vol}(ii,1:length_il) = [ seismic.trace_ilxl_bytes(1,pkey_loc) seismic.trace_ilxl_bytes(end,pkey_loc) seismic.trace_ilxl_bytes(1,skey_loc) seismic.trace_ilxl_bytes(end,skey_max_loc) seismic.trace_ilxl_bytes(1,tkey_loc) seismic.trace_ilxl_bytes(end,tkey_loc) uint64(files_in.names{il})];
            else 
                length_il = size([ seismic.trace_ilxl_bytes(1,pkey_loc) seismic.trace_ilxl_bytes(end,pkey_loc) seismic.trace_ilxl_bytes(1,skey_loc) seismic.trace_ilxl_bytes(end,skey_max_loc) uint64(files_in.names{il}) ],2);
                file_pkey_start{i_vol}(ii,1:length_il) = [ seismic.trace_ilxl_bytes(1,pkey_loc) seismic.trace_ilxl_bytes(end,pkey_loc) seismic.trace_ilxl_bytes(1,skey_loc) seismic.trace_ilxl_bytes(end,skey_max_loc) uint64(files_in.names{il}) ];
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
    job_meta.n_samples{i_vol} = seismic.n_samples;  % Number f samples in the current Volume
    job_meta.trc_head{i_vol} = 240;                 % Length of trace header??
    job_meta.bytes_per_sample{i_vol} = 4;           % Number of bytes per sample. Usually 4
    job_meta.vol_nblocks(i_vol,1) = ii-1;
    job_meta.pkey_min(i_vol,1) = min(pkey_min);     % Minimum of Primary Key (inline generally)
    job_meta.pkey_max(i_vol,1) = max(pkey_max);     % Maximum of Primary Key (inline generally)
    job_meta.pkey_inc(i_vol,1) = mode(pkey_inc);    % Mode of Primary Key (inline generally)
    job_meta.skey_min(i_vol,1) = min(skey_min);     % Minimum of Secondary Key (xline generally)
    job_meta.skey_max(i_vol,1) = max(skey_max);     % Maximum of Secondary Key (xline generally)
    job_meta.skey_inc(i_vol,1) = mode(skey_inc);    % Mode of Secondary Key (xline generally)
    
    % If Data is an angle Gather rather than angle stacks
    if seismic.is_gather == 1
        job_meta.tkey_min(i_vol,1) = min(tkey_min); % Minimum of tertiary Key (angle)
        job_meta.tkey_max(i_vol,1) = max(tkey_max); % Maximum of tertiary Key (angle)
        job_meta.tkey_inc(i_vol,1) = mode(tkey_inc);% Mode of tertiary Key (angle)
        job_meta.is_gather = 1;
    else
        job_meta.is_gather = 0;
    end  
end
%% code to find missing traces and remake mat lite files

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
% JS we should find unique pkeys for gathers files, but need to loop over block instead of vol

if is_gather == 0
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
            %cur_output_filename{i_key} = 'hfysorth';
            for file_loop = 1:1:size(file_pkey_start{1,i_vol},1)
                %if compress_ilxl_bytes(i_key,pkey_loc) > file_pkey_start{1,i_vol}(file_loop,1) && compress_ilxl_bytes(i_key,pkey_loc) < file_pkey_start{1,i_vol}(file_loop,2);
                if compress_ilxl_bytes(i_key,pkey_loc) >  min([file_pkey_start{1,i_vol}(file_loop,1) file_pkey_start{1,i_vol}(file_loop,2)]) && compress_ilxl_bytes(i_key,pkey_loc) < max([file_pkey_start{1,i_vol}(file_loop,1) file_pkey_start{1,i_vol}(file_loop,2)]);    
                    cur_output_filename{i_key} = [char(file_pkey_start{1,i_vol}(file_loop,5:end))];
                %elseif compress_ilxl_bytes(i_key,pkey_loc) == file_pkey_start{1,i_vol}(file_loop,1) || compress_ilxl_bytes(i_key,pkey_loc) == file_pkey_start{1,i_vol}(file_loop,2);
                elseif compress_ilxl_bytes(i_key,pkey_loc) == min([file_pkey_start{1,i_vol}(file_loop,1) file_pkey_start{1,i_vol}(file_loop,2)]) || compress_ilxl_bytes(i_key,pkey_loc) == max([file_pkey_start{1,i_vol}(file_loop,1) file_pkey_start{1,i_vol}(file_loop,2)]);    
                    % need to test for the skey being in the range of the file
                    % or maybe the next file
                    test_cur_skey_min = compress_ilxl_bytes(i_key,skey_loc);
                    test_cur_skey_max = compress_ilxl_bytes(i_key,skey_max_loc);
                    %if file_pkey_start{1,i_vol}(file_loop,3) <= test_cur_skey_min && file_pkey_start{1,i_vol}(file_loop,4) >= test_cur_skey_max
                    if min([file_pkey_start{1,i_vol}(file_loop,3) file_pkey_start{1,i_vol}(file_loop,4)]) <= test_cur_skey_min && max([file_pkey_start{1,i_vol}(file_loop,3) file_pkey_start{1,i_vol}(file_loop,4)])  >= test_cur_skey_max    
                        %if file_pkey_start{1,i_vol}(file_loop,3) >= test_cur_skey_min && file_pkey_start{1,i_vol}(file_loop,4) <= test_cur_skey_max
                        cur_output_filename{i_key} = [char(file_pkey_start{1,i_vol}(file_loop,5:end))];
                    else
                        %if strcmp(cur_output_filename{i_key},'hfysorth');
                            cur_output_filename{i_key} = 'hfysorth';
                        %end
                    end
                else
                    %if strcmp(cur_output_filename{i_key},'hfysorth');
                        cur_output_filename{i_key} = 'hfysorth';
                    %end
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
        
        for f_key = 1:1:size(byte_loc2,1) % this might not be looping through all required files? is this size on the wrong thing! meaning 437 is correct.
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
end
% do we need to do any of this for single volumes!
% ##################################################################
 if is_gather == 0
       job_meta.files = job_meta.files'; 
%     job_meta.files = regexprep(job_meta.files,'mat_orig_lite','mat_lite')';
end
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

dispstrj = 'starting to make blocks';
disp(dispstrj);

[job_meta.block_keys,job_meta.n_blocks] = segy_make_blocks(job_meta_path);

save(job_meta_path,'-struct','job_meta','-v7.3'); % Saves Seismic structure to mat file

% ##################################################################
% Find live blocks

job_meta.liveblocks = select_live_blocks(job_meta_path);
save(job_meta_path,'-struct','job_meta','-v7.3'); % Saves Seismic structure to mat file

fprintf('Saved seismic structure to file ...\n')

cd(start_point)
end