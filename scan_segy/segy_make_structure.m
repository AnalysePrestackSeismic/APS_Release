function segy_make_structure(filepath,il_byte,xl_byte,offset_byte,anggath,filename)
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
%   segy_make_strcuture: makes _mat_lite and mat_orig_lite files from the Seismic Header 
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
start_point = pwd;                                          % Store startpoint as goto flag to reference later on

%-----------CHECK FOR AN EXISTING SCAN---------------------------------
cd(filepath)                                                % Change directory into given file path
file_mat = regexp(filename,'\.','split');                   % Split the name into 'file name' and 'extension'
file_mat1 = strcat(filepath,file_mat{1},'.mat_orig_lite');  % Append the file name with new extension '.mat_orig_lite'

scan_exist1 = exist(file_mat1,'file');                      % Check whether this file already exists in the directory

file_mat2 = strcat(filepath,file_mat{1},'.mat_lite');       % Append the file name with new extension '.mat_lite'
scan_exist2 = exist(file_mat2,'file');                      % Check whether this file already exists in the directory

cd(start_point)                                             % Change directory to start point goto flag
usefix = 0;
%------------- CONVERT STRING (ARGUMENTS) TO DOUBLE----------------------

il_byte = str2double(il_byte);                                  % In-line byte location from argument
xl_byte = str2double(xl_byte);                                  % X-line byte location from argument
offset_byte = str2double(offset_byte);                          % Offset byte location from argument
anggath = str2double(anggath);                                  % 1 is an angle gather, other wise offset and it just indexes from 1 taken from argument

%------------ SCAN THE SEGY FILE HEADERS AND REFORMATTING------------------
if scan_exist1 == 0 && scan_exist2 == 0                         %  If no previous scan exists
    file_mat = strcat(filepath,file_mat{1},'.mat_orig_lite');   %  Append the file name with new extension '.mat_orig_lite' ??
    % Scan segy as previous scan has not been made
    % Read input file
    seismic.filepath = strcat(filepath,filename);               % Store file path in  a structure 'seismic'
    seismic.fid=fopen(seismic.filepath,'r','b');                % Find and store file ID after openning in Big-endian ordering in  a structure 'seismic'
    
    fprintf('Scan does not exist\n');                           % Tell user that scan does not exist
    fprintf('Scanning file: %s ...\n',filename);                % Tell user that you are scanning file
    
    % Read textual header (3200 bytes)
    seismic.text_header = fread(seismic.fid,3200,'uchar');      % Read the seismic EBCDIC HEADER first 3200 characters
    seismic.text_header = char(ebcdic2ascii(reshape(seismic.text_header,80,40)')); % Reshape the string of 3200 characters into a 80 x 40 matrix and convert to ASCII
    
    % Read binary header as
    seismic.binary_header = fread(seismic.fid,400,'uint8');     % Read the Binary Header
    
    % Re-interpret binary header as uint16 or uint32 as required
    two_bytes = seismic.binary_header(1:2:399)*256 + seismic.binary_header(2:2:400);
    four_bytes = ((seismic.binary_header(1:4:9)*256 + seismic.binary_header(2:4:10))*256+seismic.binary_header(3:4:11))*256+seismic.binary_header(4:4:12);
    seismic.binary_header = [four_bytes(1:3);two_bytes(7:200)];
    
    %-------------------------------------------------------------------
    %-----WRITE SOME MORE INFORMATION IN THE STUCTURE FROM HEADERS------------
    
    seismic.n_samples = seismic.binary_header(8);               % Number of Samples
    seismic.s_rate = seismic.binary_header(6);                  % Sampling Rate
    seismic.file_type = seismic.binary_header(10);              % File Type
    
    if seismic.file_type < 1 || seismic.file_type > 5 % need to break if not 1 or 5 because we don't handle it
        seismic.file_type = input('Non standard seismic file type. Please enter seismic file type (1 (IBM),2,3,4 or 5 (IEEE)): ', 's'); % Request user to manually enter file type
        seismic.file_type = str2num(seismic.file_type);         % COnvert File Type to a number
    else
        bytes_per_sample = 4;                                   % Assign default 4 bytes per sample
    end
    
    seismic.ilxl_bytes = [il_byte xl_byte];                     % Inline and X- Line Byte Locations
    %offset_byte = 37;
    % --Calculate number of traces--
    ll = dir(seismic.filepath);
    seismic.n_traces = (1/bytes_per_sample)*(ll.bytes-3600)/(seismic.n_samples+60); % Number of traces = (file size - header size)/ (size of each trace+size of trace header) 
    %-------------------------------------------------------------------
    %--------------------INTITIATE BLOCK DIVISION---------------------
    
    % Set number of traces per block
    blocktr = 10000;                                            % Edit this number if you want to change size of blocks          
    if seismic.n_traces < blocktr
        blocktr = seismic.n_traces;
    end
    loop_end = floor(seismic.n_traces/blocktr);                 % Number of blocks (Use to intiat for loops later on )
    %--------------------------------------------------------------------
    
    % -----------START WRITTING INTO .MAT_ORIG_LITE------------------------
    
    filepath_binary = uint64(seismic.filepath);                 % Convert File path into 64-bit integer array
    pad_filepath = zeros(1,(2000-length(filepath_binary)));     % Padding for file path by zeros at the end to make its length 2000
    filepath_binary = [filepath_binary,pad_filepath];           % Pad File path
    fid_write = fopen(file_mat,'w');                            % Open .mat_orig_lite for writting
    fwrite(fid_write,[filepath_binary';seismic.file_type;seismic.s_rate;seismic.n_samples;seismic.n_traces;il_byte;xl_byte;offset_byte],'double'); % write basic information into .mat_orig_lite file
    %---------------------------------------------------------------------
    %--------------READ SEGY DATA AND WRITE INFO IN .MAT_ORIG_LITE--------    
    skip_textual_binary = 3600;                                 % Length of EBCIDIC header (bytes)
    trc_head = 240;                                             % Length of Trace Header (bytes)
    trc_length = seismic.n_samples*bytes_per_sample;            % Length of trace (bytes)
    last_byte = 0;
    
    
    % Loop to read segy data 
    tic                                                         % Start clock
    for ii = 1:loop_end
        % Read blocktr x trace headers and trace data as uint32 into a temporary matrix
        tmptr = fread(seismic.fid,[120+(bytes_per_sample/2)*seismic.n_samples,blocktr],'uint16=>uint16');
        tmptr = tmptr(1:120,:);
        [trace_header bytes_to_samples] = interpret(tmptr);
        
        % Inline / crossline compression step
        trace_ilxl_bytes(:,1) = trace_header(bytes_to_samples == il_byte,:)';
        trace_ilxl_bytes(:,2) = trace_header(bytes_to_samples == xl_byte,:)';
        
        % check to see if there is a zero in the xline header
        % if so loop through all the traces and move the header location by
        % one cj99
%         if sum(trace_ilxl_bytes(:,2) == 0) > 0
%              %trace_ilxl_bytes(1,1); 
%              usefix = 1;
%              duffindxs = find(trace_ilxl_bytes(:,2) == 0);
%              %trace_ilxl_bytes(duffindxs,[1 2 4]);
%              trace_ilxl_bytes(duffindxs,1) = trace_header(bytes_to_samples == 185,duffindxs)';
%              trace_ilxl_bytes(duffindxs,2) = trace_header(bytes_to_samples == 189,duffindxs)';
%              %trace_ilxl_bytes(duffindxs,[1 2 4]);
%         end
        %
        % a check to find a particular error at an inline xline loc
        %         if sum(trace_ilxl_bytes(:,1) == 8263) > 0
        %             if sum(trace_ilxl_bytes(:,2) == 8225) > 0
        %                 trace_ilxl_bytes(1,2)
        %             end
        %         end
               
        trace_ilxl_bytes(:,3) = last_byte+(trc_head:trc_head+trc_length:blocktr*(trc_length+trc_head));
        last_byte = trace_ilxl_bytes(end,3)+trc_length; % store to add on during loop
        trace_ilxl_bytes(:,3) = trace_ilxl_bytes(:,3)+skip_textual_binary;
        
        % Check for gathers
        n_traces_to_check = 1000;
        if n_traces_to_check > blocktr
            n_traces_to_check = blocktr;
        end
        % test to see if there are any duplicate inline xline locations in
        % the first 1000 traces, if yes then a gather if no then angle
        % stacks
        if length(unique(trace_ilxl_bytes(1:n_traces_to_check,1:2),'rows')) < n_traces_to_check
            % is gathers
            is_gather = 1;
            if ii == 1
                fwrite(fid_write,is_gather,'double');
            end
            trace_ilxl_bytes(:,4) = trace_header(bytes_to_samples == offset_byte,:)'; % offset hard wired
            %                 if ii == 32
%             if usefix == 1 % cj99
%                 trace_ilxl_bytes(duffindxs,4) =  trace_header(bytes_to_samples == 35,duffindxs)';
%                 usefix = 0;
%             end
            %                 end
            if anggath == 1
                compress_ilxl_bytes = gather_compress_ilxl_bytes(trace_ilxl_bytes,blocktr);
            else
                compress_ilxl_bytes = gather_compress_ilxl_bytes_offset(trace_ilxl_bytes,blocktr);
            end
            % compress_ilxl_bytes(1:(end-1),9) = (compress_ilxl_bytes(2:end,3)-compress_ilxl_bytes(1:(end-1),3)) ./  ( (((compress_ilxl_bytes(1:(end-1),7)- compress_ilxl_bytes(1:(end-1),6))./compress_ilxl_bytes(1:(end-1),8))+1) .*((seismic.n_samples*4)+240));
        else
            is_gather = 0;
            if ii == 1
                fwrite(fid_write,is_gather,'double');
            end
            
            compress_ilxl_bytes = trace_compress_ilxl_bytes(trace_ilxl_bytes,blocktr);
        end
        
        fwrite(fid_write,reshape(compress_ilxl_bytes',[],1),'double');
        usefix = 0;
        % fprintf('Block %d completed...\n',ii);
    end
    
    % Calculate the number of trace headers not read by the loop above
    leftovers = seismic.n_traces-loop_end*blocktr;
    
    % Read the remaining trace headers (if any)
    clearvars tmptrheader trace_header bytes_to_samples trace_ilxl_bytes
    if leftovers > 0
        tmptr = fread(seismic.fid,[120+2*seismic.n_samples,leftovers],'uint16=>uint16');
        tmptr = tmptr(1:120,:);
        [trace_header bytes_to_samples] = interpret(tmptr);
        % tmptrheader(1:120,1+loop_end*blocktr:seismic.n_traces) = tmptr(1:120,:);
        
        trace_ilxl_bytes(:,1) = trace_header(bytes_to_samples == il_byte,:)';
        trace_ilxl_bytes(:,2) = trace_header(bytes_to_samples == xl_byte,:)';
        trace_ilxl_bytes(:,3) = last_byte+(trc_head:trc_head+trc_length:leftovers*(trc_length+trc_head));
        trace_ilxl_bytes(:,3) = trace_ilxl_bytes(:,3)+skip_textual_binary;
        
%         % fix for a particular error in a segy file cj99
%         if sum(trace_ilxl_bytes(:,2) == 0) > 0
%              %trace_ilxl_bytes(1,1); 
%              usefix = 1;
%              duffindxs = find(trace_ilxl_bytes(:,2) == 0);
%              trace_ilxl_bytes(duffindxs,1) = trace_header(bytes_to_samples == 185,duffindxs)';
%              trace_ilxl_bytes(duffindxs,2) = trace_header(bytes_to_samples == 189,duffindxs)';
%         end
        
        % Check for gathers
        if is_gather == 1;
            % is gathers
            trace_ilxl_bytes(:,4) = trace_header(bytes_to_samples == offset_byte,:)'; % offset hard wired
            
%             if usefix == 1 % cj99
%                 trace_ilxl_bytes(duffindxs,4) =  trace_header(bytes_to_samples == 35,duffindxs)';
%                 usefix = 0;
%             end
            
            if anggath == 1
                compress_ilxl_bytes = gather_compress_ilxl_bytes(trace_ilxl_bytes,leftovers);
            else
                compress_ilxl_bytes = gather_compress_ilxl_bytes_offset(trace_ilxl_bytes,leftovers);
            end
        else
            compress_ilxl_bytes = trace_compress_ilxl_bytes(trace_ilxl_bytes,leftovers);
        end
        
        fwrite(fid_write,reshape(compress_ilxl_bytes',[],1),'double');
        usefix = 0;
    end
    
    %         for i_live = 1:1:size(non_live_traces,2)
    %             non_live_traces{i_live} = non_live_traces{i_live}+(blocktr*i_live-1);
    %         end
    %         non_live_traces = cell2mat(non_live_traces');
    
    fprintf('%d MB of data read (and written) at %d MB/sec in %d seconds\n',round((ll.bytes-3600)/(1024*1024)),round(((ll.bytes-3600)/(1024*1024))/toc),round(toc));
    
    % Close segy file
    
    fclose('all');
    
else
    
    fprintf('\nScan of segy %s already exists\n',filename);
    
end
clearvars tmptrheader trace_header bytes_to_samples trace_ilxl_bytes
end


function [trace_header bytes_to_samples] = interpret(tmptrheader)
byte_type = [ ...
    2*ones(7,1); ones(4,1);
    2*ones(8,1); ones(2,1);
    2*ones(4,1); ones(46,1);
    2*ones(5,1); ones(2,1);
    2*ones(1,1); ones(5,1);
    2*ones(1,1); ones(1,1);
    2*ones(1,1); ones(2,1);
    2*ones(1,1); 2*ones(1,1)];

ntr = size(tmptrheader,2);
trace_header = zeros(91,ntr);
bytes_to_samples = zeros(91,1);

count =1;
for ii = 1:91
    bytes_to_samples(ii,1) = 2*count-1;
    if byte_type(ii) == 1
        trace_header(ii,:) = double(tmptrheader(count,:));
        count = count+1;
    elseif byte_type(ii) == 2
        trace_header(ii,:) = double(tmptrheader(count,:))*2^16 + double(tmptrheader(count+1,:));
        count = count+2;
    end
end

trace_header(21,:) = trace_header(21,:)-2^16;

end

function ascii=ebcdic2ascii(ebcdic)
% Function converts EBCDIC string to ASCII
% see http://www.room42.com/store/computer_center/code_tables.shtml
%
% Written by: E. Rietsch: Feb. 20, 2000
% Last updated:
%
%           ascii=ebcdic2ascii(ebcdic)
% INPUT
% ebcdic    EBCDIC string
% OUTPUT
% ascii	   ASCII string

pointer= ...
    [ 0    16    32    46    32    38    45    46    46    46    46    46   123   125    92    48
    1    17    33    46    46    46    47    46    97   106   126    46    65    74    46    49
    2    18    34    50    46    46    46    46    98   107   115    46    66    75    83    50
    3    19    35    51    46    46    46    46    99   108   116    46    67    76    84    51
    4    20    36    52    46    46    46    46   100   109   117    46    68    77    85    52
    5    21    37    53    46    46    46    46   101   110   118    46    69    78    86    53
    6    22    38    54    46    46    46    46   102   111   119    46    70    79    87    54
    7    23    39    55    46    46    46    46   103   112   120    46    71    80    88    55
    8    24    40    56    46    46    46    46   104   113   121    46    72    81    89    56
    9    25    41    57    46    46    46    46   105   114   122    46    73    82    90    57
    10    26    42    58    46    33   124    58    46    46    46    46    46    46    46    46
    11    27    43    59    46    36    44    35    46    46    46    46    46    46    46    46
    12    28    44    60    60    42    37    64    46    46    46    46    46    46    46    46
    13    29    45    61    40    41    95    39    46    46    91    93    46    46    46    46
    14    30    46    46    43    59    62    61    46    46    46    46    46    46    46    46
    15    31    47    63   124    94    63    34    46    46    46    46    46    46    46    46];

pointer=reshape(pointer,1,256);

ascii=pointer(ebcdic+1);

end

% function compress_ilxl_bytes  = trace_compress_ilxl_bytes(trace_ilxl_bytes,blocktr)
%
% pkey_loc = 1; % column numbers needs to be implemented
% skey_loc = 2;
% byte_loc = 3;
% tkey_loc = 4;
% tkey_max_loc = 5;
% tkey_inc_loc = 6;
%
% start_idx = 1;
% count = 0;
% row_i = 1;
%
% %blocktr = size(trace_ilxl_bytes,1);
% pkey_prev = -995837;
% skey_prev = -9999437;
% skey_inc = -999971;
% cur_inc = -27389995;
%
% %compress_ilxl_bytes = zeros([blocktr,5],'int64');
%
% if blocktr > 1
%     for row_i = start_idx:blocktr
%         pkey = trace_ilxl_bytes(row_i,pkey_loc);
%         skey = trace_ilxl_bytes(row_i,skey_loc);
%         tbyte = trace_ilxl_bytes(row_i,byte_loc);
%
%         if pkey == pkey_prev
%             cur_inc = skey - skey_prev;
%             if cur_inc ~= skey_inc
%                 if cur_inc == 0
%                     count = count + 1;
%                     compress_ilxl_bytes(count,:) = [ pkey skey tbyte skey 1 ];
%                     skey_inc = -999971;
%                     skey_prev = skey;
%                 else
%                     if skey_inc == -999971 % cur_inc is first time or after duplicate trace
%                         compress_ilxl_bytes(count,4:5) = [ skey cur_inc ];
%                         skey_inc = cur_inc;
%                         skey_prev = skey;
%                     else % cur_inc is not 0
%                         count = count + 1;
%                         compress_ilxl_bytes(count,:) = [ pkey skey tbyte skey cur_inc ];
%                         skey_inc = cur_inc;
%                         skey_prev = skey;
%                     end
%                 end
%             else % cur_inc == skey_inc
%                 compress_ilxl_bytes(count,4) = skey;
%                 skey_prev = skey;
%             end
%
%
%         else % pkey ~= pkey_prev
%             count = count + 1;
%             compress_ilxl_bytes(count,:) = [ pkey skey tbyte skey 1 ];
%             pkey_prev = pkey;
%             skey_prev = skey;
%             skey_inc = -999971;
%         end
%
%     end
%
% else % for blocktr = 1
%     count = count + 1;
%     compress_ilxl_bytes(count,:) = [ trace_ilxl_bytes(row_i,pkey_loc) trace_ilxl_bytes(row_i,skey_loc) trace_ilxl_bytes(row_i,byte_loc) trace_ilxl_bytes(row_i,skey_loc) 1 ];
% end
% % now truncate the array
% %compress_ilxl_bytes(count+1:end,:) = [];
%
% end

% function compress_ilxl_bytes = gather_compress_ilxl_bytes(trace_ilxl_bytes,blocktr)
%     pkey_loc = 1; % column numbers needs to be implemented
%     skey_loc = 2;
%     byte_loc = 3;
%     tkey_loc = 4;
%
%     trace_ilxl_bytes(1:end-1,5) = diff(trace_ilxl_bytes(:,4));
%
%     start_idx = 1;
%     count = 1;
%     row_i = 1;
%     pkey_prev = -995837;
%     skey_prev = -9999437;
%
%     if blocktr > 1
%         for row_i = start_idx:blocktr
%             pkey = trace_ilxl_bytes(row_i,pkey_loc);
%             skey = trace_ilxl_bytes(row_i,skey_loc);
%             tkey = trace_ilxl_bytes(row_i,tkey_loc);
%             tbyte = trace_ilxl_bytes(row_i,byte_loc);
%             tkey_inc = trace_ilxl_bytes(row_i,5);
%
%             if pkey == pkey_prev % same inline
%                 if tkey_inc == tkey_inc_prev
%                     tkey_inc_prev = tkey_inc;
%                     skey_prev = skey;
%                 elseif skey == skey_prev
%                     compress_offset_bytes(count,5:6) = [tkey tkey_inc_prev];
%                     count = count + 1;
%                     tkey_inc_prev = tkey_inc;
%                     skey_prev = skey;
%                 else
%                     compress_offset_bytes(count,:) = [ pkey skey tbyte tkey tkey 1];
%                     tkey_inc_prev = tkey_inc;
%                     skey_prev = skey;
%                 end
%
%             else % pkey ~= pkey_prev
%                 %count = count + 1;
%                 compress_offset_bytes(count,:) = [ pkey skey tbyte tkey tkey 1];
%                 pkey_prev = pkey;
%                 skey_prev = skey;
%                 tkey_inc_prev = tkey_inc;
%             end
%
%         end
%
%     % Now compress inlines, crosslines with same offset range
%     blocktr = size(compress_offset_bytes,1);
%
%     if blocktr > 1
%         start_idx = 1;
%         count = 0;
%         row_i = 1;
%
%         %blocktr = size(trace_ilxl_bytes,1);
%         pkey_prev = -995837;
%         skey_prev = -9999437;
%         skey_inc = -999971;
%         cur_inc = -27389995;
%         tkey_min_prev = -999971;
%         tkey_max_prev = 999971;
%         tkey_inc_prev = 999971;
%         tkey_inc_prev = -27389995;
%
%         %compress_ilxl_bytes = zeros([blocktr,5],'int64');
%             for row_i = start_idx:blocktr
%                 pkey = compress_offset_bytes(row_i,pkey_loc);
%                 skey = compress_offset_bytes(row_i,skey_loc);
%                 tbyte = compress_offset_bytes(row_i,byte_loc);
%                 tkey_min = compress_offset_bytes(row_i,4);
%                 tkey_max = compress_offset_bytes(row_i,5);
%                 tkey_inc = compress_offset_bytes(row_i,6);
%
%                 if pkey == pkey_prev && tkey_min == tkey_min_prev && tkey_max == tkey_max_prev && tkey_inc == tkey_inc_prev
%                     cur_inc = skey - skey_prev;
%                     if cur_inc ~= skey_inc
%                         if cur_inc == 0
%                             count = count + 1;
%                             compress_ilxl_bytes(count,:) = [ pkey skey tbyte skey 1 tkey_min tkey_max tkey_inc];
%                             skey_inc = -999971;
%                             skey_prev = skey;
%                             tkey_min_prev = tkey_min;
%                             tkey_max_prev = tkey_max;
%                             tkey_inc_prev = tkey_inc;
%                         else
%                             if skey_inc == -999971 % cur_inc is first time or after duplicate trace
%                                 compress_ilxl_bytes(count,4:5) = [ skey cur_inc ];
%                                 skey_inc = cur_inc;
%                                 skey_prev = skey;
%                                 tkey_min_prev = tkey_min;
%                                 tkey_max_prev = tkey_max;
%                                 tkey_inc_prev = tkey_inc;
%                             else % cur_inc is not 0
%                                 count = count + 1;
%                                 compress_ilxl_bytes(count,:) = [ pkey skey tbyte skey cur_inc tkey_min tkey_max tkey_inc];
%                                 skey_inc = cur_inc;
%                                 skey_prev = skey;
%                                 tkey_min_prev = tkey_min;
%                                 tkey_max_prev = tkey_max;
%                                 tkey_inc_prev = tkey_inc;
%                             end
%                         end
%                     else % cur_inc == skey_inc
%                         compress_ilxl_bytes(count,4) = skey;
%                         skey_prev = skey;
%                         tkey_min_prev = tkey_min;
%                         tkey_max_prev = tkey_max;
%                         tkey_inc_prev = tkey_inc;
%                     end
%
%
%                 else % pkey ~= pkey_prev
%                     count = count + 1;
%                     compress_ilxl_bytes(count,:) = [ pkey skey tbyte skey 1 tkey_min tkey_max tkey_inc];
%                     pkey_prev = pkey;
%                     skey_prev = skey;
%                     tkey_min_prev = tkey_min;
%                     tkey_max_prev = tkey_max;
%                     tkey_inc_prev = tkey_inc;
%                     skey_inc = -999971;
%                 end
%
%             end
%     else
%         count = count + 1;
%         compress_ilxl_bytes(count,:) = [ compress_offset_bytes(row_i,pkey_loc) compress_offset_bytes(row_i,skey_loc) compress_offset_bytes(row_i,byte_loc) compress_offset_bytes(row_i,skey_loc) 1 compress_offset_bytes(row_i,4) compress_offset_bytes(row_i,5) compress_offset_bytes(row_i,6)];
%     end
%
%     else % for blocktr = 1
%         count = count + 1;
%         compress_ilxl_bytes(count,:) = [ trace_ilxl_bytes(row_i,pkey_loc) trace_ilxl_bytes(row_i,skey_loc) trace_ilxl_bytes(row_i,byte_loc) trace_ilxl_bytes(row_i,skey_loc) 1 trace_ilxl_bytes(row_i,4) trace_ilxl_bytes(row_i,4) 1];
%     end
%
% end

