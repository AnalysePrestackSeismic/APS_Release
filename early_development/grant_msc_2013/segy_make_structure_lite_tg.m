function [seismic] = segy_make_structure_lite_tg(filepath,filename,il_byte,xl_byte,varargin)
%%-------------------------------------------------------------------------
%   Run this on one angle stack as they should all have the same
%   geometry/file structure. Typical IL/XL byte locations: 189 & 193. This
%   will create a mat file in the location of the input SEGY definining the
%   structure of the input file.
%%-------------------------------------------------------------------------

    start_point = pwd; % remember starting directory
    
    seismic.filepath = strcat(filepath,filename);
    
    % Check for existing scan
    cd(filepath)
    
    file_mat = regexp(filename,'\.','split');
    file_mat = strcat(file_mat{1},'.mat');
    scan_exist = exist(file_mat,'file');
    
    cd(start_point)
       
    if scan_exist == 0
        % Scan segy as previous scan has not been made
        % Read input file
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
        
        % Calculate number of traces
        ll = dir(seismic.filepath);
        seismic.n_traces = 0.25*(ll.bytes-3600)/(seismic.n_samples+60);
        
        % Set number of traces per block
        blocktr = 2000;
        loop_end = floor(seismic.n_traces/blocktr);
        
        % Loop to read segy data
        tic
        for ii = 1:loop_end
            % Read blocktr x trace headers and trace data as uint32 into a temporary matrix
            tmptr = fread(seismic.fid,[120+2*seismic.n_samples,blocktr],'uint16=>uint16');
            tmptrheader(1:120,1+(ii-1)*blocktr:ii*blocktr) = tmptr(1:120,:);
        end
        
        % Calculate the number of trace headers not read by the loop above
        leftovers = seismic.n_traces-loop_end*blocktr;

        % Read the remaining trace headers (if any)
        if leftovers > 0
            tmptr = fread(seismic.fid,[120+2*seismic.n_samples,leftovers],'uint16=>uint16');
            tmptrheader(1:120,1+loop_end*blocktr:seismic.n_traces) = tmptr(1:120,:);
        end
        fprintf('%d MB of data read at %d MB/sec in %d seconds\n',round((ll.bytes-3600)/(1024*1024)),round(((ll.bytes-3600)/(1024*1024))/toc),round(toc));

        [seismic.trace_header seismic.bytes_to_samples] = interpret(tmptrheader);
        fprintf('Interpreting headers ...\n')
        
        % Extract inline - crossline information
        seismic.ilxl_bytes = [il_byte xl_byte]; 
        seismic.trace_ilxl_bytes(:,1) = seismic.trace_header(seismic.bytes_to_samples == il_byte,:)';
        seismic.trace_ilxl_bytes(:,2) = seismic.trace_header(seismic.bytes_to_samples == xl_byte,:)';

        coscaler = seismic.trace_header(21,1);
        if coscaler < 0
            coscaler = 1/abs(coscaler);
        elseif coscaler == 0
            coscaler = 1;
        end
        
        % extra bytes 
        seismic.extra_bytes_data = 0;       
        if (size(varargin,2) > 0)
            seismic.extra_bytes = cell2mat(varargin);
            seismic.extra_bytes_data = zeros(seismic.n_traces,length(seismic.extra_bytes));            
            for ii=1:1:length(seismic.extra_bytes) 
               seismic.extra_bytes_data(:,ii) = seismic.trace_header(seismic.bytes_to_samples == seismic.extra_bytes(ii),:)';
               if or(and(seismic.extra_bytes(ii)>=73,seismic.extra_bytes(ii)<=85),and(seismic.extra_bytes(ii)>=181,seismic.extra_bytes(ii)<=185))
                   seismic.extra_bytes_data(:,ii) = seismic.extra_bytes_data(:,ii)*coscaler;
               end
            end            
        end
        
        % Create trace byte pointers
        skip_textual_binary = 3600;
        trc_head = 240;
        trc_length = seismic.n_samples*4;
        seismic.trace_ilxl_bytes(1:1:seismic.n_traces,3) = ....
                 skip_textual_binary +...
                 cumsum(repmat(trc_head,1,seismic.n_traces)) +... 
                 (0:1:seismic.n_traces-1).*trc_length;
             
        c = clock;
        seismic.file_id = c*c';  
        
        seismic.trc_head_length = trc_head;
       
        % to be able to exclude traces outside of range
        seismic.min_iline = min(seismic.trace_ilxl_bytes(:,1));
        seismic.max_iline = max(seismic.trace_ilxl_bytes(:,1));
        seismic.min_xline = min(seismic.trace_ilxl_bytes(:,2));
        seismic.max_xline = max(seismic.trace_ilxl_bytes(:,2));

        % find inline/crossline increments
        seismic.n_iline = length(unique(seismic.trace_ilxl_bytes(:,1)));
        seismic.n_xline = length(unique(seismic.trace_ilxl_bytes(:,2)));
        seismic.il_inc = floor(mean(diff(unique(seismic.trace_ilxl_bytes(:,1)))));
        seismic.xl_inc = floor(mean(diff(unique(seismic.trace_ilxl_bytes(:,2)))));
        
        % calculate number of crosslines per inline
        ils = unique(seismic.trace_ilxl_bytes(:,1));
        loop = 1;
        xl_per_il = zeros(length(ils),1);
        for ii=1:seismic.n_traces
            if ils(loop,1) == seismic.trace_ilxl_bytes(ii,1);
                xl_per_il(loop,1) = xl_per_il(loop,1)+1;
            else
                loop=loop+1;
                xl_per_il(loop,1) = 1;
            end
        end
        
        % calculate survey boundary inline/crossline/byte positions
        last_xl_per_il = cumsum(xl_per_il);
        first_xl_per_il = last_xl_per_il-xl_per_il+1;       
        for ii=1:length(ils)
            seismic.boundary(ii,:) = [ils(ii,1),seismic.trace_ilxl_bytes(first_xl_per_il(ii),2:3),ils(ii,1),seismic.trace_ilxl_bytes(last_xl_per_il(ii),2:3)];
        end        
        seismic.boundary = reshape(seismic.boundary',3,[])';
        
        seismic.il_start_bytes = [seismic.boundary(logical(repmat([1;0],length(ils),1)),:),xl_per_il];
        
        % Close segy file
        fclose(seismic.fid);
        
        % Scan scan to /mat file
        cd(filepath)
        fprintf('Saving seismic structure to file ...\n')
        save(file_mat,'-struct','seismic','-v7.3')
        fwrite(fopen(strcat(file_mat,'_lite'),'w'),[seismic.n_samples;reshape(seismic.trace_ilxl_bytes',[],1)],'double');
        cd(start_point)
        fprintf('Scan of file saved as %s\n',file_mat);
        fclose all;
    else
        % Scan of segy already exists so open it
        fprintf('\nScan of segy %s already exists\n',filename);
        fprintf('Opening .mat file %s\n',file_mat);
        cd(filepath)
        seismic = load(file_mat,'-mat');
        fwrite(fopen(strcat(file_mat,'_lite'),'w'),[seismic.n_samples;reshape(seismic.trace_ilxl_bytes',[],1)],'double');
        cd(start_point)        
    end 
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