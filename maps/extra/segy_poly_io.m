function segy_poly_io(xloc_byte,yloc_byte,filename)
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
% segy_poly_io(xloc_byte,yloc_byte,filename) 
% this reads a xy polygon and zeros out all the data that is not inside
% that polygon and then writes it back out to a new ieee segy file
%%-------------------------------------------------------------------------
% constants

% Tanzania Block 3 Polygon
xv = [638079.01
616591.06
606476.20
641453.01
672954.15
683680.74
668103.64
649928.29
644000.99
639741.26
636100.63
638079.01];

yv = [9021703.90
9016716.11
9060287.78
9060171.21
9060066.22
9013682.83
9009988.54
9023908.23
9049263.22
9030813.33
9029965.06
9021703.90];

% xv = [446915.03
%     449919.10
%     449873.95
%     451878.93
%     451835.28
%     451769.89
%     451726.31
%     451682.74
%     451639.19
%     449624.21
%     449646.89
%     448639.86
%     448663.00
%     448686.14
%     447680.02
%     447703.62
%     447727.23
%     446722.02
%     446770.15
%     444761.55
%     444786.53
%     443782.69
%     443808.13
%     446818.31
%     446915.03];

% yv = [ 6362402.15
%     6362364.24
%     6358652.98
%     6358628.90
%     6354918.14
%     6349352.10
%     6345641.43
%     6341930.78
%     6338220.14
%     6338244.28
%     6340099.59
%     6340112.03
%     6341967.34
%     6343822.65
%     6343835.33
%     6345690.64
%     6347545.96
%     6347558.88
%     6351269.53
%     6351296.08
%     6353151.40
%     6353165.04
%     6355020.37
%     6354980.19
%     6362402.15];

xv = xv.*100;
yv = yv.*100;
xv = int32(xv);
yv = int32(yv);

blocktr = 20000;                                        % Edit this number if you want to change size of blocks
usefix = 0;
%------------- CONVERT STRING (ARGUMENTS) TO DOUBLE----------------------

xloc_byte = str2double(xloc_byte);                                  % In-line byte location from argument
yloc_byte = str2double(yloc_byte);                                  % X-line byte location from argument
%coscal = str2double(coscal);                          % Offset byte location from argument
%scale = str2double(scale);                                      % scale the data with an envelope scaling

%------------ SCAN THE SEGY FILE HEADERS AND REFORMATTING------------------


% Read input file
strseismic.filepath = filename;               % Store file path in  a structure 'strseismic'
strseismic.fid=fopen(strseismic.filepath,'r','b');                % Find and store file ID after openning in Big-endian ordering in  a structure 'strseismic'

%fprintf('Scan does not exist\n');                           % Tell user that scan does not exist
%fprintf('Scanning file: %s ...\n',filename);                % Tell user that you are scanning file

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

out_file_name = strcat(filename,'_subpoly_ieee.segy');
out_file_name_slice = strcat(filename,'_slice');
%out_file_name = fullfile(out_dir, sprintf('%s_block_%d%s',results_out{i_results,1},i_block,'.segy'));

%% Initialise files and headers
poly_ieee_segy = fopen(out_file_name,'w'); % Open file for writing
% Write out the ebcdic_header
fwrite(poly_ieee_segy,strseismic.text_header_orig,'uchar',0,'ieee-be');

strseismic.binary_header_orig(26) = 5;

% Write out the binary header
fwrite(poly_ieee_segy,strseismic.binary_header_orig,'uint8',0,'ieee-be');

%-----WRITE SOME MORE INFORMATION IN THE STUCTURE FROM HEADERS------------

strseismic.n_samples = strseismic.binary_header(8);               % Number of Samples
strseismic.s_rate = strseismic.binary_header(6);                  % Sampling Rate
strseismic.file_type = strseismic.binary_header(10);              % File Type

if strseismic.file_type < 1 || strseismic.file_type > 5 % need to break if not 1 or 5 because we don't handle it
    strseismic.file_type = input('Non standard strseismic file type. Please enter strseismic file type (1 (IBM),2,3,4 or 5 (IEEE)): ', 's'); % Request user to manually enter file type
    strseismic.file_type = str2num(strseismic.file_type);         % COnvert File Type to a number
else
    bytes_per_sample = 4;                                   % Assign default 4 bytes per sample
end

strseismic.ilxl_bytes = [xloc_byte yloc_byte];                     % Inline and X- Line Byte Locations
%offset_byte = 37;
% --Calculate number of traces--
ll = dir(strseismic.filepath);
strseismic.n_traces = floor((1/bytes_per_sample)*(ll.bytes-3600)/(strseismic.n_samples+60)); % Number of traces = (file size - header size)/ (size of each trace+size of trace header)
int_n_traces = int32(strseismic.n_traces);
%--------------------INTITIATE BLOCK DIVISION---------------------

% Set number of traces per block

if strseismic.n_traces < blocktr
    blocktr = strseismic.n_traces;
end
loop_end = floor(strseismic.n_traces/blocktr);                 % Number of blocks (Use to intiat for loops later on )

qcslice = zeros(6,int_n_traces,'single');               % decalre the array for a couple of time slices

%--------------------------------------------------------------------


%--------------READ SEGY DATA --------
% skip_textual_binary = 3600;                                 % Length of EBCIDIC header (bytes)
% trc_head = 240;                                             % Length of Trace Header (bytes)
% trc_length = strseismic.n_samples*bytes_per_sample;            % Length of trace (bytes)
% last_byte = 0;
% dhexa = uint32(hex2dec('7f000000'));
% dhexb = uint32(hex2dec('00ffffff'));
% cjconst = double(1/(2^24));
% cj16arr = double(ones(strseismic.n_samples,blocktr).*16);

% Loop to read segy data
tic                                                         % Start clock
for ii = 1:loop_end
    % Read blocktr x trace headers and trace data as uint32 into a temporary matrix
    if strseismic.file_type == 1
        % Convert traces from IBM32FP read as UINT32 into IEEE64FP (doubles) - need to make it singles
        traces_tmp = fread(strseismic.fid,[60+strseismic.n_samples,blocktr],'*uint32');
        
        trchead = traces_tmp(1:60,:);
        [trace_header bytes_to_samples] = interpretbe(reshape(typecast(trchead(:),'uint16'),120,[]));
        
        traces = single((1-2*double(bitget(traces_tmp(61:end,:),32))).*16.^ ...
            (double(bitshift(bitand(traces_tmp(61:end,:),uint32(hex2dec('7f000000'))),-24))-64).* ...
            (double(bitand(traces_tmp(61:end,:),uint32(hex2dec('00ffffff'))))/2^24));
        
%         % tried lots of profiling but nothing was quicker then the compound
%         % funcion, perhaps compiled versions would be
%         tracesb = single((1-2*double(bitget(traces_tmp(61:end,:),32))).*16.^(double(bitshift(bitand(traces_tmp(61:end,:),dhexa),-24))-64).*(double(bitand(traces_tmp(61:end,:),dhexb))*cjconst));
%         %tracesb2 = single((1-2*double(bitget(traces_tmp(61:end,:),32))).*realpow(16,(double(bitshift(bitand(traces_tmp(61:end,:),dhexa),-24))-64)).*(double(bitand(traces_tmp(61:end,:),dhexb))*cjconst));
%         %tracesb2 = single((1-2*double(bitget(traces_tmp(61:end,:),32))).*bsxfun(@power,cj16arr,(double(bitshift(bitand(traces_tmp(61:end,:),dhexa),-24))-64)).*(double(bitand(traces_tmp(61:end,:),dhexb))*cjconst));
%         
%         %         tracesc1 = 1-2*double(bitget(traces_tmp(61:end,:),32));
%         %         %tracesc2 = 16.^(double(bitshift(bitand(traces_tmp(61:end,:),dhexa),-24))-64);
%         %         %tracesc2a = 16;
%         %         tracesc2b = (double(bitshift(bitand(traces_tmp(61:end,:),dhexa),-24))-64);
%         %         %tracesc2 = 16.^tracesc2b;
%         %         %tracesc2 = realpow(cj16arr,tracesc2b);
%         %         tracesc2 = bsxfun(@power,cj16arr,tracesc2b);
%         %
%         %         tracesc3 = double(bitand(traces_tmp(61:end,:),dhexb))*cjconst;
%         %
%         %         tracesc = single(tracesc1.*tracesc2.*tracesc3);
        
    elseif strseismic.file_type == 2
        disp('This strseismic file type is not currently supported. Please speak to Charles Jones.');
    elseif strseismic.file_type == 5
        % Traces are IEEE32FP (singles)
        traces = fread(strseismic.fid,[60+strseismic.n_samples,blocktr],strcat(num2str(strseismic.n_samples),'*float32=>float32'));
        %trace_headers = typecast(single(reshape(traces(1:60,:),1,60*n_traces_to_read)),'int32');
        %trace_headers = reshape(trace_headers,60,n_traces_to_read);
        
        trchead = typecast(single(reshape(traces(1:60,:),1,60*blocktr)),'uint16');
        %trchead = reshape(trchead,120,n_traces_to_read);
        
        [trace_header bytes_to_samples] = interpretbe(reshape(trchead,120,[]));
        
        traces = traces(61:end,:);
    else
        disp('This strseismic file type is not currently supported. Please speak to Charles Jones.');
    end
    
    % Inline / crossline reading step
    trace_ilxl_bytes(:,1) = int32(trace_header(bytes_to_samples == xloc_byte,:))';
    trace_ilxl_bytes(:,2) = int32(trace_header(bytes_to_samples == yloc_byte,:))';
    %offset_read = int32(trace_header(bytes_to_samples == strseismic.tkey,:))';
    
    % make a couple of matlab arrays to hold a couple of time slices
    qcslice(1:4,(((ii-1)*blocktr)+1):ii*blocktr) = [single(trace_ilxl_bytes(:,1)');single(trace_ilxl_bytes(:,2)');traces(floor((60+strseismic.n_samples)/3),:);traces(floor((60+strseismic.n_samples)/2),:)];
    
    %close;
    %-------------- zero out the traces not in the polygon --------------
    %calculate if trace headers are in polygon or not
    tr_use = inpolygon(trace_ilxl_bytes(:,1),trace_ilxl_bytes(:,2),xv,yv);
    %use logical to set the third dimension of the trace length to zero
    tracesout = bsxfun(@times,traces,tr_use');
    
    if sum(tr_use) > 0
        % make a couple of matlab arrays to hold a couple of time slices
        qcslice(5:6,(((ii-1)*blocktr)+1):ii*blocktr) = [tracesout(floor((60+strseismic.n_samples)/3),:);tracesout(floor((60+strseismic.n_samples)/2),:)];
        
        %reattach the raw headers with a reshape to float32
        temparr2 = [reshape(typecast(trchead(:),'int32'),size(trchead,1),blocktr); reshape(typecast(tracesout(:),'int32'),strseismic.n_samples,blocktr)];
        
        %and write out as big endian ieee segy, have to write ebcdic first though
        fwrite(poly_ieee_segy,temparr2,'int32',0,'ieee-be');
    end
   
    
end
clearvars trace_header bytes_to_samples trace_ilxl_bytes tr_use tracesout trchead temparr2

%-------------------------- read the remaining traces --------------------
leftovers = strseismic.n_traces-(loop_end*blocktr);

if leftovers > 0
    %         tmptr = fread(strseismic.fid,[120+2*strseismic.n_samples,leftovers],'uint16=>uint16');
    % Read blocktr x trace headers and trace data as uint32 into a temporary matrix
    if strseismic.file_type == 1
        clearvars traces_tmp traces 
        % Convert traces from IBM32FP read as UINT32 into IEEE64FP (doubles) - need to make it singles
        traces_tmp = fread(strseismic.fid,[60+strseismic.n_samples,leftovers],'*uint32');
        
        trchead = traces_tmp(1:60,:);
        [trace_header bytes_to_samples] = interpretbe(reshape(typecast(trchead(:),'uint16'),120,[]));
        
        traces = single((1-2*double(bitget(traces_tmp(61:end,:),32))).*16.^ ...
            (double(bitshift(bitand(traces_tmp(61:end,:),uint32(hex2dec('7f000000'))),-24))-64).* ...
            (double(bitand(traces_tmp(61:end,:),uint32(hex2dec('00ffffff'))))/2^24));
        
    elseif strseismic.file_type == 2
        disp('This strseismic file type is not currently supported. Please speak to Charles Jones.');
    elseif strseismic.file_type == 5
        % Traces are IEEE32FP (singles)
        clearvars traces 
        traces = fread(strseismic.fid,[60+strseismic.n_samples,leftovers],strcat(num2str(strseismic.n_samples),'*float32=>float32'));
        %trace_headers = typecast(single(reshape(traces(1:60,:),1,60*n_traces_to_read)),'int32');
        %trace_headers = reshape(trace_headers,60,n_traces_to_read);
        
        trchead = typecast(single(reshape(traces(1:60,:),1,60*leftovers)),'uint16');
        %trchead = reshape(trchead,120,n_traces_to_read);
        
        [trace_header bytes_to_samples] = interpretbe(reshape(trchead,120,[]));
        
        traces = traces(61:end,:);
    else
        disp('This strseismic file type is not currently supported. Please speak to Charles Jones.');
    end
    
    % Inline / crossline reading step
    trace_ilxl_bytes(:,1) = int32(trace_header(bytes_to_samples == xloc_byte,:))';
    trace_ilxl_bytes(:,2) = int32(trace_header(bytes_to_samples == yloc_byte,:))';
    %
    % make a couple of matlab arrays to hold a couple of time slices
    qcslice(1:4,(((ii)*blocktr)+1):end) = [single(trace_ilxl_bytes(:,1)');single(trace_ilxl_bytes(:,2)');traces(floor((60+strseismic.n_samples)/3),:);traces(floor((60+strseismic.n_samples)/2),:)];
    
    %close;
    %-------------- zero out the traces not in the polygon --------------
    tr_use = inpolygon(trace_ilxl_bytes(:,1),trace_ilxl_bytes(:,2),xv,yv);
    tracesout = bsxfun(@times,traces,tr_use');
    
    if sum(tr_use) > 0
        % make a couple of matlab arrays to hold a couple of time slices
        qcslice(5:6,(((ii)*blocktr)+1):end) = [tracesout(floor((60+strseismic.n_samples)/3),:);tracesout(floor((60+strseismic.n_samples)/2),:)];
        
        %reattach the raw headers with a reshape to float32
        temparr2 = [reshape(typecast(trchead(:),'int32'),size(trchead,1),leftovers); reshape(typecast(tracesout(:),'int32'),strseismic.n_samples,leftovers)];
        
        %and write out as big endian ieee segy, have to write ebcdic first though
        fwrite(poly_ieee_segy,temparr2,'int32',0,'ieee-be');
    end
end

fclose(strseismic.fid);
fclose(poly_ieee_segy);
toc
fprintf('%d MB of data read (and some written) at %d MB/sec in %d seconds\n',round((ll.bytes-3600)/(1024*1024)),round(((ll.bytes-3600)/(1024*1024))/toc),round(toc));


% save data as a mat file and then display plots using another routine
endcj = strseismic.n_traces;
out_file_name_slice_tmp = strcat(out_file_name_slice,'_qc_out.mat');
save(out_file_name_slice_tmp,'qcslice','endcj','-v7.3'); % Saves Seismic structure to mat file

fprintf('wrote out file %s \n now view using plot_segy_poly_qc_slices(''%s'')\n',out_file_name_slice_tmp,out_file_name_slice_tmp);


% red = [repmat(1,50,1);(1:-0.02:0)'];
% blue = [(0:0.02:1)';repmat(1,50,1)];
% green = [(0:0.02:1)';(0.98:-0.02:0)'];
% redblue = [red green blue];
% 
% scrsz = get(0,'ScreenSize');
% cjf = figure('Position',[1 scrsz(4) scrsz(3)/2 scrsz(4)],'PaperType','A1','Visible','off','PaperPositionMode','auto','PaperOrientation','landscape');
% colormap(redblue);
% endcj = strseismic.n_traces;
% colrngcj = (0.5*max(qcslice(3,1:2:endcj)));
% %props = {'LineStyle','none','Marker','o','MarkerEdge','b','MarkerSize',6};
% %line([qcslice(1,1:10:endcj),qcslice(1,1:10:endcj)],[qcslice(2,1:10:endcj),qcslice(2,1:10:endcj)],props{:});
% scatter(qcslice(1,1:10:endcj),qcslice(2,1:10:endcj),20,qcslice(3,1:10:endcj),'filled');
% caxis([(-1*colrngcj) colrngcj]);
% 
% 
% out_file_name_slice_tmp = strcat(out_file_name_slice,'_whole.pdf');
% 
% print(cjf,'-dpdf','-noui','-r200',out_file_name_slice_tmp);
% 
% 
% scatter(qcslice(1,1:10:endcj),qcslice(2,1:10:endcj),20,qcslice(5,1:10:endcj),'filled');
% caxis([(-1*colrngcj) colrngcj]);
% out_file_name_slice_tmp = strcat(out_file_name_slice,'_sub_polygon.pdf');
% print(cjf,'-dpdf','-noui','-r200',out_file_name_slice_tmp);

%close cjf;






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
