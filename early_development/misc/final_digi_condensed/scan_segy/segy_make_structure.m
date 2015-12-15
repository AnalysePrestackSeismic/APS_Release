function segy_make_structure(filepath,il_byte,xl_byte,filename)
%function [seismic] = segy_make_structure(filepath,filename,il_byte,xl_byte)
%%-------------------------------------------------------------------------
%   Run this on one angle stack as they should all have the same
%   geometry/file structure. Typical IL/XL byte locations: 189 & 193. This
%   will create a mat file in the location of the input SEGY definining the
%   structure of the input file.
%
%   Structure format
%   PKey SKey Byte_Loc SKey_max SKey_inc TKey TKey_max TKey_inc
%
%   Use segy_index_checker(seismic_mat_path) to check scan validity
%%-------------------------------------------------------------------------
    start_point = pwd;
    
    % Check for existing scan
    cd(filepath)
    file_mat = regexp(filename,'\.','split');
    file_mat = strcat(filepath,file_mat{1},'.mat_orig_lite');
    scan_exist = exist(file_mat,'file');

    cd(start_point)
    
    il_byte = str2double(il_byte);
    xl_byte = str2double(xl_byte);

    if scan_exist == 0
        % Scan segy as previous scan has not been made
        % Read input file
        seismic.filepath = strcat(filepath,filename);
        seismic.fid=fopen(seismic.filepath,'r','b');

        fprintf('Scan does not exist\n');
        fprintf('Scanning file: %s ...\n',filename);

        % Read textual header (3200 bytes)
        seismic.text_header = fread(seismic.fid,3200,'uchar');
        seismic.text_header = char(ebcdic2ascii(reshape(seismic.text_header,80,40)'));    

        % Read binary header as 
        seismic.binary_header=fread(seismic.fid,400,'uint8');

        % Re-interpret binary header as uint16 or uint32 as required
        two_bytes=seismic.binary_header(1:2:399)*256+seismic.binary_header(2:2:400);
        four_bytes=((seismic.binary_header(1:4:9)*256+seismic.binary_header(2:4:10))*256+seismic.binary_header(3:4:11))*256+seismic.binary_header(4:4:12);
        seismic.binary_header=[four_bytes(1:3);two_bytes(7:200)];

        seismic.n_samples = seismic.binary_header(8);
        seismic.s_rate = seismic.binary_header(6);
        seismic.file_type = seismic.binary_header(10);

        if seismic.file_type < 1 || seismic.file_type > 5 % need to break if not 1 or 5 because we don't handle it
            seismic.file_type = input('Non standard seismic file type. Please enter seismic file type (1 (IBM),2,3,4 or 5 (IEEE)): ', 's');
            seismic.file_type = str2num(seismic.file_type);
        else
            bytes_per_sample = 4;
        end

        seismic.ilxl_bytes = [il_byte xl_byte]; 
        offset_byte = 37;
        % Calculate number of traces
        ll = dir(seismic.filepath);
        seismic.n_traces = (1/bytes_per_sample)*(ll.bytes-3600)/(seismic.n_samples+60);

        % Set number of traces per block
        blocktr = 10000;
        loop_end = floor(seismic.n_traces/blocktr);

        filepath_binary = uint64(seismic.filepath);
        pad_filepath = zeros(1,(2000-length(filepath_binary)));
        filepath_binary = [filepath_binary,pad_filepath];
        fid_write = fopen(file_mat,'w');        
        fwrite(fid_write,[filepath_binary';seismic.file_type;seismic.s_rate;seismic.n_samples;seismic.n_traces;il_byte;xl_byte;offset_byte],'double');  

        skip_textual_binary = 3600;
        trc_head = 240;
        trc_length = seismic.n_samples*bytes_per_sample;
        last_byte = 0; 
        
        % Loop to read segy data        
        tic
        for ii = 1:loop_end
            % Read blocktr x trace headers and trace data as uint32 into a temporary matrix
            tmptr = fread(seismic.fid,[120+(bytes_per_sample/2)*seismic.n_samples,blocktr],'uint16=>uint16');
            tmptr = tmptr(1:120,:);
            [trace_header bytes_to_samples] = interpret(tmptr);            

            % Inline / crossline compression step
            trace_ilxl_bytes(:,1) = trace_header(bytes_to_samples == il_byte,:)';
            trace_ilxl_bytes(:,2) = trace_header(bytes_to_samples == xl_byte,:)';

            trace_ilxl_bytes(:,3) = last_byte+(trc_head:trc_head+trc_length:blocktr*(trc_length+trc_head));
            last_byte = trace_ilxl_bytes(end,3)+trc_length; % store to add on during loop                         
            trace_ilxl_bytes(:,3) = trace_ilxl_bytes(:,3)+skip_textual_binary;

            % Check for gathers
            if length(unique(trace_ilxl_bytes(1:1000,1:2),'rows')) < 1000 
                % is gathers
                is_gather = 1;
                if ii == 1
                    fwrite(fid_write,is_gather,'double'); 
                end

                trace_ilxl_bytes(:,4) = trace_header(bytes_to_samples == offset_byte,:)'; % offset hard wired
                compress_ilxl_bytes = gather_compress_ilxl_bytes(trace_ilxl_bytes,blocktr);
            else
                is_gather = 0;
                if ii == 1
                    fwrite(fid_write,is_gather,'double'); 
                end

                compress_ilxl_bytes = trace_compress_ilxl_bytes(trace_ilxl_bytes,blocktr);
            end                                               

            fwrite(fid_write,reshape(compress_ilxl_bytes',[],1),'double');      

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
            % Check for gathers
            if is_gather == 1;
                % is gathers
                trace_ilxl_bytes(:,4) = trace_header(bytes_to_samples == offset_byte,:)'; % offset hard wired
                compress_ilxl_bytes = gather_compress_ilxl_bytes(trace_ilxl_bytes,leftovers);
            else
                compress_ilxl_bytes = trace_compress_ilxl_bytes(trace_ilxl_bytes,leftovers);
            end 

            fwrite(fid_write,reshape(compress_ilxl_bytes',[],1),'double');

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

function [compress_ilxl_bytes]  = trace_compress_ilxl_bytes_orig(trace_ilxl_bytes,blocktr)
    pkey_loc = 1; % column numbers needs to be implemented
    skey_loc = 2;
    byte_loc = 3;
    skey_max_loc = 4;
    skey_inc_loc = 5;
    tkey_loc = 6;
    tkey_max_loc = 7;
    tkey_inc_loc = 8;

    start_idx = 1;
    count = 1;
    constant = 10000;
    
%     % Calculate live traces
%     pkey_min = min(trace_ilxl_bytes(:,pkey_loc));
%     pkey_max = max(trace_ilxl_bytes(:,pkey_loc));
%     pkey_inc = mode(diff(unique(trace_ilxl_bytes(:,pkey_loc))));
%     if isnan(pkey_inc)
%         pkey_inc = 1;
%         iln = 1;
%     else
%         iln = 1+((pkey_max - pkey_min)/pkey_inc);
%     end
%     skey_min = min(trace_ilxl_bytes(:,skey_loc));
%     skey_max = max(trace_ilxl_bytes(:,skey_loc));
%     skey_inc = mode(diff(unique(trace_ilxl_bytes(:,skey_loc))));     
%     if isnan(skey_inc)
%         skey_inc = 1;
%         xln = 1;
%     else    
%         xln = 1+((skey_max - skey_min)/skey_inc);
%     end
% 
%     bounding_box = [reshape(repmat((pkey_min:pkey_inc:pkey_max),xln,1),[],1),...
%     repmat((skey_min:skey_inc:skey_max)',iln,1),(1:1:iln*xln)'];
% 
%     % clculate live traces as a loop
%     % Get the linear indices of the lives traces in the survey bounding box
%     live_traces = bounding_box(ismember(bounding_box(:,1:2),trace_ilxl_bytes(:,1:2),'rows'),3);
%     
%     non_live_traces = find(diff(live_traces) > 1);
    
    row_i = 0;
    while start_idx < blocktr                    
            cdp_1 = (trace_ilxl_bytes(start_idx,1)+constant)*(trace_ilxl_bytes(start_idx,2));
        if start_idx < blocktr-1
            cdp_2 = (trace_ilxl_bytes(start_idx+1,1)+constant)*(trace_ilxl_bytes(start_idx+1,2));
        end
        if start_idx < blocktr-2
            cdp_3 = (trace_ilxl_bytes(start_idx+2,1)+constant)*(trace_ilxl_bytes(start_idx+2,2));
        end
        if (cdp_2 - cdp_1) == (cdp_3 - cdp_2)                          
            row_i = row_i+1;            
        elseif trace_ilxl_bytes(start_idx,1) == trace_ilxl_bytes(start_idx+1,1)  

            compress_ilxl_bytes(count,1) = trace_ilxl_bytes(start_idx-row_i,1);
            compress_ilxl_bytes(count,2) = trace_ilxl_bytes(start_idx-row_i,2);
            compress_ilxl_bytes(count,3) = trace_ilxl_bytes(start_idx-row_i,3);             
            compress_ilxl_bytes(count,4) = trace_ilxl_bytes(start_idx+1,2);

            xl_inc = trace_ilxl_bytes(start_idx-row_i+1,2)-trace_ilxl_bytes(start_idx-row_i,2);
            compress_ilxl_bytes(count,5) = xl_inc;           

            row_i = 0;
            count = count + 1;
        end
        start_idx = start_idx + 1;
    end
    compress_ilxl_bytes(end,4) = trace_ilxl_bytes(end,2);   
end

function [compress_ilxl_bytes]  = trace_compress_ilxl_bytes(trace_ilxl_bytes,blocktr)

pkey_loc = 1; % column numbers needs to be implemented
skey_loc = 2;
byte_loc = 3;
skey_max_loc = 4;
skey_inc_loc = 5;
tkey_loc = 6;
tkey_max_loc = 7;
tkey_inc_loc = 8;

start_idx = 1;
count = 0;
row_i = 1;

%blocktr = size(trace_ilxl_bytes,1);
pkey_prev = -995837;
skey_prev = -9999437;
skey_inc = -999971;
cur_inc = -27389995;

%compress_ilxl_bytes = zeros([blocktr,5],'int64');

if blocktr > 1
    for row_i = start_idx:blocktr
        pkey = trace_ilxl_bytes(row_i,pkey_loc);
        skey = trace_ilxl_bytes(row_i,skey_loc);
        tbyte = trace_ilxl_bytes(row_i,byte_loc);
        
        if pkey == pkey_prev
            cur_inc = skey - skey_prev;
            if cur_inc ~= skey_inc
                if cur_inc == 0
                    count = count + 1;
                    compress_ilxl_bytes(count,:) = [ pkey skey tbyte skey 1 ];
                    skey_inc = -999971;
                    skey_prev = skey;
                else
                    if skey_inc == -999971 % cur_inc is first time or after duplicate trace
                        compress_ilxl_bytes(count,4:5) = [ skey cur_inc ];
                        skey_inc = cur_inc;
                        skey_prev = skey;
                    else % cur_inc is not 0
                        count = count + 1;
                        compress_ilxl_bytes(count,:) = [ pkey skey tbyte skey cur_inc ];
                        skey_inc = cur_inc;
                        skey_prev = skey;
                    end
                end
            else % cur_inc == skey_inc
                compress_ilxl_bytes(count,4) = skey;
                skey_prev = skey;
            end
            
            
        else % pkey ~= pkey_prev
            count = count + 1;
            compress_ilxl_bytes(count,:) = [ pkey skey tbyte skey 1 ];
            pkey_prev = pkey;
            skey_prev = skey;
            skey_inc = -999971;
        end
        
    end
    
else % for blocktr = 1
    count = count + 1;
    compress_ilxl_bytes(count,:) = [ trace_ilxl_bytes(row_i,pkey_loc) trace_ilxl_bytes(row_i,skey_loc) trace_ilxl_bytes(row_i,byte_loc) trace_ilxl_bytes(row_i,skey_loc) 1 ];
end
% now truncate the array
%compress_ilxl_bytes(count+1:end,:) = [];

end



function compress_ilxl_bytes = gather_compress_ilxl_bytes(trace_ilxl_bytes,blocktr)
    start_idx = 1;
    count = 1;   
    row_i = 0;
    while start_idx < blocktr                  
            offset_1 = trace_ilxl_bytes(start_idx,4);
        if start_idx < blocktr-1            
            offset_2 = trace_ilxl_bytes(start_idx+1,4);
        end
        if start_idx < blocktr-2            
            offset_3 = trace_ilxl_bytes(start_idx+2,4);
        end
        if (offset_2 - offset_1) == (offset_3-offset_2)                          
            row_i = row_i+1;            
        elseif offset_3 < offset_2   
            compress_offset_bytes(count,1) = trace_ilxl_bytes(start_idx-row_i,1);
            compress_offset_bytes(count,2) = trace_ilxl_bytes(start_idx-row_i,2);
            compress_offset_bytes(count,3) = trace_ilxl_bytes(start_idx-row_i,3);           

            compress_offset_bytes(count,4) = trace_ilxl_bytes(start_idx-row_i,4);
            compress_offset_bytes(count,5) = trace_ilxl_bytes(start_idx+1,4);
            
            offset_inc = trace_ilxl_bytes(start_idx-row_i+1,4)-trace_ilxl_bytes(start_idx-row_i,4);
            compress_offset_bytes(count,6) = offset_inc;

            row_i = 0;
            count = count + 1;
        end
        start_idx = start_idx + 1;
    end 
    start_idx = start_idx - 1;
    compress_offset_bytes(count,1) = trace_ilxl_bytes(start_idx-row_i,1);
    compress_offset_bytes(count,2) = trace_ilxl_bytes(start_idx-row_i,2);
    compress_offset_bytes(count,3) = trace_ilxl_bytes(start_idx-row_i,3);
    compress_offset_bytes(count,4) = trace_ilxl_bytes(start_idx-row_i,4);
    compress_offset_bytes(count,5) = trace_ilxl_bytes(start_idx+1,4);
    offset_inc = trace_ilxl_bytes(start_idx-row_i+1,4)-trace_ilxl_bytes(start_idx-row_i,4);
    compress_offset_bytes(count,6) = offset_inc;
    
    % check for inline crossline patterns
    start_idx = 1;
    count = 1;  
    row_i = 0;
    constant = 10000;
    n_traces = size(compress_offset_bytes,1);
    while start_idx < n_traces            
            cdp_1 = (compress_offset_bytes(start_idx,1)+constant)*(compress_offset_bytes(start_idx,2));
            soffset_1 = compress_offset_bytes(start_idx,4);
            eoffset_1 = compress_offset_bytes(start_idx,5);
        if start_idx < n_traces-1
            cdp_2 = (compress_offset_bytes(start_idx+1,1)+constant)*(compress_offset_bytes(start_idx+1,2));
            soffset_2 = compress_offset_bytes(start_idx+1,4);
            eoffset_2 = compress_offset_bytes(start_idx+1,5);
        end
        if start_idx < n_traces-2
            cdp_3 = (compress_offset_bytes(start_idx+2,1)+constant)*(compress_offset_bytes(start_idx+2,2));
            soffset_3 = compress_offset_bytes(start_idx+2,4);
            eoffset_3 = compress_offset_bytes(start_idx+2,5);
        end
        if (cdp_2 - cdp_1) == (cdp_3 - cdp_2) && (soffset_2 - soffset_1) == (soffset_3 - soffset_2) && (eoffset_2 - eoffset_1) == (eoffset_3 - eoffset_2)
            row_i = row_i+1;
        elseif compress_offset_bytes(start_idx,1) == compress_offset_bytes(start_idx+1,1)  
            compress_ilxl_bytes(count,1) = compress_offset_bytes(start_idx-row_i,1);
            compress_ilxl_bytes(count,2) = compress_offset_bytes(start_idx-row_i,2);
            compress_ilxl_bytes(count,3) = compress_offset_bytes(start_idx-row_i,3);             
            compress_ilxl_bytes(count,4) = compress_offset_bytes(start_idx+1,2);

            xl_inc = compress_offset_bytes(start_idx-row_i+1,2)-compress_offset_bytes(start_idx-row_i,2);
            compress_ilxl_bytes(count,5) = xl_inc; 
            
            compress_ilxl_bytes(count,6) = compress_offset_bytes(start_idx-row_i,4);
            compress_ilxl_bytes(count,7) = compress_offset_bytes(start_idx+1,5);           
            
            compress_ilxl_bytes(count,8) = compress_offset_bytes(start_idx-row_i,6);

            row_i = 0;
            count = count + 1;            
        end
        start_idx = start_idx + 1;
    end
    if start_idx == n_traces
        compress_ilxl_bytes(count,1) = compress_offset_bytes(end,1);
        compress_ilxl_bytes(count,2) = compress_offset_bytes(end,2);
        compress_ilxl_bytes(count,3) = compress_offset_bytes(end,3);
        compress_ilxl_bytes(count,4) = compress_offset_bytes(end,2);
        compress_ilxl_bytes(count,5) = compress_ilxl_bytes(count,4)-compress_ilxl_bytes(count,2);
        compress_ilxl_bytes(count,6) = compress_offset_bytes(end,4);
        compress_ilxl_bytes(count,7) = compress_offset_bytes(end,5);
        compress_ilxl_bytes(count,8) = compress_offset_bytes(end,6);
    end

end

