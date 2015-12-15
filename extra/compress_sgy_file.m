function compress_sgy_file(filename,compress_type)
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
% compress_segy
% this reads segy file and compresses it and writes it out
% to a new ieee segy file with scalars attached to the bottom of the
% traces
%%-------------------------------------------------------------------------
% constants
blocktr = 500;
%blocktr = 20000;                                        % Edit this number if you want to change size of blocks
usefix = 0;
%------------- CONVERT STRING (ARGUMENTS) TO DOUBLE----------------------

compress_type = str2double(compress_type);                                  % compression type argument
% xloc_byte = str2double(xloc_byte);                                  % In-line byte location from argument
% yloc_byte = str2double(yloc_byte);                                  % X-line byte location from argument
%------------ SCAN THE SEGY FILE HEADERS AND REFORMATTING------------------

% Read input file
strseismic.filepath = filename;               % Store file path in  a structure 'strseismic'
strseismic.fid=fopen(strseismic.filepath,'r','b');                % Find and store file ID after openning in Big-endian ordering in  a structure 'strseismic'

fprintf('reading file: %s \n',filename);                % Tell user that you are scanning file

% Read textual header (3200 bytes)
strseismic.text_header_orig = fread(strseismic.fid,3200,'uchar');      % Read the strseismic EBCDIC HEADER first 3200 characters
strseismic.text_header = strseismic.text_header_orig;
strseismic.text_header = char(ebcdic2ascii(reshape(strseismic.text_header,80,40)')); % Reshape the string of 3200 characters into a 80 x 40 matrix and convert to ASCII

% Read binary header as
strseismic.binary_header_orig = fread(strseismic.fid,400,'uint8');     % Read the Binary Header
strseismic.binary_header = strseismic.binary_header_orig;

% Re-interpret binary header as uint16 or uint32 as required
two_bytes = strseismic.binary_header(1:2:399)*256 + strseismic.binary_header(2:2:400);
four_bytes = ((strseismic.binary_header(1:4:9)*256 + strseismic.binary_header(2:4:10))*256+strseismic.binary_header(3:4:11))*256+strseismic.binary_header(4:4:12);
strseismic.binary_header = [four_bytes(1:3);two_bytes(7:200)];


%------------------------- start writting the output file ------------------------------------------

out_file_name = strcat(filename,'_comp.segy');
fprintf('outputting to file: %s \n',out_file_name); 
%out_file_name_slice = strcat(filename,'_slice');

% Initialise files and headers -----------------------------------------
comp_ieee_sgy = fopen(out_file_name,'w'); % Open file for writing
% Write out the ebcdic_header
fwrite(comp_ieee_sgy,strseismic.text_header_orig,'uchar',0,'ieee-be');

strseismic.n_samples = strseismic.binary_header(8);               % Number of Samples

% set the binary to be either 16 or 8 bit
if compress_type < 3
    strseismic.binary_header_orig(26) = 8;
    extrazeros =  4 - ((strseismic.n_samples/4 - floor(strseismic.n_samples/4)) * 4);
    output_samples = strseismic.n_samples+extrazeros; 
    strseismic.binary_header_orig(21) = floor(output_samples*0.00390625);
    strseismic.binary_header_orig(22) =  (output_samples*0.00390625 - floor(output_samples*0.00390625) )* 256;
else
    strseismic.binary_header_orig(26) = 3;
    extrazeros =  2 - ((strseismic.n_samples/2 - floor(strseismic.n_samples/2)) * 2);
    output_samples = strseismic.n_samples+extrazeros; 
    strseismic.binary_header_orig(8) = output_samples;    
end

% Write out the binary header
fwrite(comp_ieee_sgy,strseismic.binary_header_orig,'uint8',0,'ieee-be');

%-----WRITE SOME MORE INFORMATION IN THE STUCTURE FROM HEADERS------------


strseismic.s_rate = strseismic.binary_header(6);                  % Sampling Rate
strseismic.file_type = strseismic.binary_header(10);              % File Type

if strseismic.file_type < 1 || strseismic.file_type > 5 % need to break if not 1 or 5 because we don't handle it
    strseismic.file_type = input('Non standard strseismic file type. Please enter strseismic file type (1 (IBM),2,3,4 or 5 (IEEE)): ', 's'); % Request user to manually enter file type
    strseismic.file_type = str2num(strseismic.file_type);         % COnvert File Type to a number
else
    bytes_per_sample = 4;                                   % Assign default 4 bytes per sample
end

% strseismic.ilxl_bytes = [xloc_byte yloc_byte];                     % Inline and X- Line Byte Locations

% --Calculate number of traces--
ll = dir(strseismic.filepath);
strseismic.n_traces = floor((1/bytes_per_sample)*(ll.bytes-3600)/(strseismic.n_samples+60)); % Number of traces = (file size - header size)/ (size of each trace+size of trace header)
int_n_traces = int64(strseismic.n_traces);
%--------------------INTITIATE BLOCK DIVISION---------------------

% Set number of traces per block

if strseismic.n_traces < blocktr
    blocktr = strseismic.n_traces;
end
loop_end = floor(strseismic.n_traces/blocktr);                 % Number of blocks (Use to intiat for loops later on )

% qcslice = zeros(6,int_n_traces,'single');               % decalre the array for a couple of time slices

%--------------READ SEGY DATA --------
% Loop to read segy data
tic                                                         % Start clock
for ii = 1:loop_end
    % Read blocktr x trace headers and trace data as uint32 into a temporary matrix
    if strseismic.file_type == 1
        % Convert traces from IBM32FP read as UINT32 into IEEE64FP (doubles) - need to make it singles
        traces_tmp = fread(strseismic.fid,[60+strseismic.n_samples,blocktr],'*uint32');
        
        trchead = traces_tmp(1:60,:);
        [trace_header, bytes_to_samples] = interpretbe(reshape(typecast(trchead(:),'uint16'),120,[]));
        
        traces = single((1-2*double(bitget(traces_tmp(61:end,:),32))).*16.^ ...
            (double(bitshift(bitand(traces_tmp(61:end,:),uint32(hex2dec('7f000000'))),-24))-64).* ...
            (double(bitand(traces_tmp(61:end,:),uint32(hex2dec('00ffffff'))))/2^24));
               
    elseif strseismic.file_type == 2
        disp('This strseismic file type is not currently supported. Please speak to Charles Jones.');
    elseif strseismic.file_type == 5
        % Traces are IEEE32FP (singles)
        traces = fread(strseismic.fid,[60+strseismic.n_samples,blocktr],strcat(num2str(strseismic.n_samples),'*float32=>float32'));       
        trchead = typecast(single(reshape(traces(1:60,:),1,60*blocktr)),'uint16');       
        [trace_header, bytes_to_samples] = interpretbe(reshape(trchead,120,[]));
        traces = traces(61:end,:);
    else
        disp('This strseismic file type is not currently supported. Please speak to Charles Jones.');
    end
    
    % Inline / crossline reading step
    %trace_ilxl_bytes(:,1) = int64(trace_header(bytes_to_samples == xloc_byte,:))';
    %trace_ilxl_bytes(:,2) = int64(trace_header(bytes_to_samples == yloc_byte,:))';
    %offset_read = int32(trace_header(bytes_to_samples == strseismic.tkey,:))';
    
    % make a couple of matlab arrays to hold a couple of time slices
    %qcslice(1:4,(((ii-1)*blocktr)+1):ii*blocktr) = [single(trace_ilxl_bytes(:,1)');single(trace_ilxl_bytes(:,2)');traces(floor((60+strseismic.n_samples)/3),:);traces(floor((60+strseismic.n_samples)/2),:)];
    
    %-------------- compress the traces --------------
    % 
    [ trim_data_filt, scalars_out, zlocs_out, tpgrad_orig, scalepeak,  orig_nsamples, percenterr ] = compress_seis_to_16_or_8bit_int_nz( traces, num2str(compress_type), '0' );
    
    % make a couple of matlab arrays to hold a couple of time slices
    %qcslice(5:6,(((ii-1)*blocktr)+1):ii*blocktr) = [tracesout(floor((60+strseismic.n_samples)/3),:);tracesout(floor((60+strseismic.n_samples)/2),:)];
    
    %reattach the raw headers with a reshape to float32
    %extrazeros =  4 - ((strseismic.n_samples/4 - floor(strseismic.n_samples/4)) * 4);
    if compress_type < 3
        trim_data_filt = [trim_data_filt; zeros(extrazeros,size(trim_data_filt,2),'int8')];
        
        temparr2 = [reshape(typecast(trchead(:),'int32'),[],blocktr);reshape(swapbytes(typecast(trim_data_filt(:),'int32')),((output_samples)/4),blocktr)];
       
        %and write out as big endian ieee segy, have to write ebcdic first though
        fwrite(comp_ieee_sgy,temparr2,'int32',0,'ieee-be');
    else
        trim_data_filt = [trim_data_filt; zeros(extrazeros,size(trim_data_filt,2),'int16')];
    end
    
    fprintf('%d traces read\n',blocktr*ii);
    
end
clearvars trace_header bytes_to_samples trace_ilxl_bytes tr_use tracesout trchead temparr2

%-------------------------- read the remaining traces --------------------
leftovers = strseismic.n_traces-(loop_end*blocktr);

if leftovers > 0
    % Read blocktr x trace headers and trace data as uint32 into a temporary matrix
    if strseismic.file_type == 1
        clearvars traces_tmp traces 
        % Convert traces from IBM32FP read as UINT32 into IEEE64FP (doubles) - need to make it singles
        traces_tmp = fread(strseismic.fid,[60+strseismic.n_samples,leftovers],'*uint32');
        
        trchead = traces_tmp(1:60,:);
        [trace_header, bytes_to_samples] = interpretbe(reshape(typecast(trchead(:),'uint16'),120,[]));
        
        traces = single((1-2*double(bitget(traces_tmp(61:end,:),32))).*16.^ ...
            (double(bitshift(bitand(traces_tmp(61:end,:),uint32(hex2dec('7f000000'))),-24))-64).* ...
            (double(bitand(traces_tmp(61:end,:),uint32(hex2dec('00ffffff'))))/2^24));
        
    elseif strseismic.file_type == 2
        disp('This strseismic file type is not currently supported. Please speak to Charles Jones.');
    elseif strseismic.file_type == 5
        % Traces are IEEE32FP (singles)
        clearvars traces 
        traces = fread(strseismic.fid,[60+strseismic.n_samples,leftovers],strcat(num2str(strseismic.n_samples),'*float32=>float32'));      
        trchead = typecast(single(reshape(traces(1:60,:),1,60*leftovers)),'uint16');
        [trace_header, bytes_to_samples] = interpretbe(reshape(trchead,120,[]));
        traces = traces(61:end,:);
    else
        disp('This strseismic file type is not currently supported. Please speak to Charles Jones.');
    end
    
    % Inline / crossline reading step
    %trace_ilxl_bytes(:,1) = int64(trace_header(bytes_to_samples == xloc_byte,:))';
    %trace_ilxl_bytes(:,2) = int64(trace_header(bytes_to_samples == yloc_byte,:))';
    %
    % make a couple of matlab arrays to hold a couple of time slices
    %qcslice(1:4,(((ii)*blocktr)+1):end) = [single(trace_ilxl_bytes(:,1)');single(trace_ilxl_bytes(:,2)');traces(floor((60+strseismic.n_samples)/3),:);traces(floor((60+strseismic.n_samples)/2),:)];
    

     %-------------- compress the traces --------------

    [ trim_data_filt, scalars_out, zlocs_out, tpgrad_orig, scalepeak,  orig_nsamples, percenterr ] = compress_seis_to_16_or_8bit_int_nz( traces, num2str(compress_type), '0' );
    
        % make a couple of matlab arrays to hold a couple of time slices
        %qcslice(5:6,(((ii)*blocktr)+1):end) = [tracesout(floor((60+strseismic.n_samples)/3),:);tracesout(floor((60+strseismic.n_samples)/2),:)];
    
    %reattach the raw headers with a reshape to float32
    %extrazeros =  4 - ((strseismic.n_samples/4 - floor(strseismic.n_samples/4)) * 4);
    if compress_type < 3
        trim_data_filt = [trim_data_filt; zeros(extrazeros,size(trim_data_filt,2),'int8')];
        
        temparr2 = [reshape(typecast(trchead(:),'int32'),[],leftovers);reshape(swapbytes(typecast(trim_data_filt(:),'int32')),((output_samples)/4),leftovers)];
       
        %and write out as big endian ieee segy, have to write ebcdic first though
        fwrite(comp_ieee_sgy,temparr2,'int32',0,'ieee-be');
    else
        trim_data_filt = [trim_data_filt; zeros(extrazeros,size(trim_data_filt,2),'int16')];
    end
       
end

fclose(strseismic.fid);
fclose(comp_ieee_sgy);
toc
fprintf('%d MB of data read (and some written) at %d MB/sec in %d seconds\n',round((ll.bytes-3600)/(1024*1024)),round(((ll.bytes-3600)/(1024*1024))/toc),round(toc));


% save data as a mat file and then display plots using another routine
% endcj = strseismic.n_traces;
% out_file_name_slice_tmp = strcat(out_file_name_slice,'_qc_out.mat');
% save(out_file_name_slice_tmp,'qcslice','endcj','-v7.3'); % Saves Seismic structure to mat file
% 
% fprintf('wrote out file %s \n now view using plot_segy_poly_qc_slices(''%s'')\n',out_file_name_slice_tmp,out_file_name_slice_tmp);



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


function [trace_header bytes_to_samples] = interpretbe(tmptrheader)
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
        trace_header(ii,:) = double(tmptrheader(count+1,:))*2^16 + double(tmptrheader(count,:)); % note this is big-endian and different to one in segy_make_structure
        count = count+2;
    end
end

trace_header(21,:) = trace_header(21,:)-2^16;

end
