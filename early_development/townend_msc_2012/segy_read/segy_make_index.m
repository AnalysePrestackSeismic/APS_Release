function [seismic]=segy_make_index(inputfile,il_byte,xl_byte,varargin)
% James Selvage, March 2012
%   
% Reference: S4M
% Need to make this read groups of traces
    
    seismic.filepaths = inputfile;
   
    % Read input file
    seismic.fid = open_segy_file(seismic.filepaths(:));
    
    fprintf('\nScanning segy %s\n',inputfile);
        
    % Read textual header (3200 bytes)
    seismic.text_header(:,:) = read_textual_file_header(seismic.fid);
        
    % Read binary header (400 bytes)
    seismic.binary_header(:,:) = read_binary_file_header(seismic.fid);
    
    % Need to vary these
    seismic.n_samples = seismic.binary_header(8);
    seismic.s_rate = seismic.binary_header(6);
              
    % Find out number of traces
    seismic.n_traces = get_no_of_traces(seismic.filepaths(:),seismic.n_samples);

    % Read trace headers (240 bytes) (file,byte_index,trace_index)
    seismic.trace_header = read_trace_headers(seismic.fid,seismic.n_traces,seismic.n_samples); 

    % Extract inline - crossline information
    index = ((il_byte-1)/4)+1; % 240 bytes is read into an array of length 60 
    seismic.trace_pointers(:,1) = get_trace_header(seismic.n_traces,seismic.trace_header,index);
    index = ((xl_byte-1)/4)+1; 
    seismic.trace_pointers(:,2) = get_trace_header(seismic.n_traces,seismic.trace_header,index);
    
    % enter bytes 
    if (size(varargin,2) > 0)
        seismic.extr_bytes = cell2mat(varargin);
        point_col = 4;
        for ii=1:1:length(seismic.extr_bytes)
           index = ((seismic.extr_bytes(ii)-1)/4)+1;  
           col = point_col+ii-1;
           seismic.trace_pointers(:,col) = get_trace_header(seismic.n_traces,seismic.trace_header,index);
        end
    end
    
    % Create trace byte pointers
    skip_textual_binary = 3600;
    trc_head = 240;
    trc_length = seismic.n_samples*4;
    seismic.trace_pointers(1:1:seismic.n_traces,3) = ....
             skip_textual_binary +...
             cumsum(repmat(trc_head,1,seismic.n_traces)) +...
             (0:1:seismic.n_traces-1).*trc_length;
         
    % to be able to exclude traces outside of range
    seismic.min_iline = min(seismic.trace_pointers(:,1));
    seismic.max_iline = max(seismic.trace_pointers(:,1));
    seismic.min_xline = min(seismic.trace_pointers(:,2));
    seismic.max_xline = max(seismic.trace_pointers(:,2));
    
    % find inline/crossline increments
    seismic.n_iline = length(unique(seismic.trace_pointers(:,1)));
    seismic.n_xline = length(unique(seismic.trace_pointers(:,2)));
    seismic.il_inc = mean(diff(unique(seismic.trace_pointers(:,1))));
    seismic.xl_inc = mean(diff(unique(seismic.trace_pointers(:,2))));
         
    % Remove trace headers     
    seismic = rmfield(seismic, 'trace_header');
    
    % Close segy file
    close_segy_file(seismic.fid);
    
   %file = regexprep(inputfile, '.', 'repstr')
%    save_input = input('Do you want to save scan to .mat file? [1 - Yes, 0 - No]: '); % (could also save to ascii?)
%    if (save_input == 1)        
%         save_name = input('Save name: ');
%         save(save_name) 
%    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fid = open_segy_file(filepaths)
   
    mformat='ieee-be';
    fid=fopen(filepaths,'r',mformat);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function close_segy_file(fid)
   
    fclose(fid);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function text_header=read_textual_file_header(fid)

    text_header=fread(fid,3200,'uchar');
    
    text_header=char(ebcdic2ascii(reshape(text_header,80,40)'));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function binary_header=read_binary_file_header(fid)
% Read binary header and initiate the seismic dataset

    bh=fread(fid,400,'uchar');

    two_bytes=bh(1:2:399)*256+bh(2:2:400);
    four_bytes=((bh(1:4:9)*256+bh(2:4:10))*256+bh(3:4:11))*256+bh(4:4:12);

    binary_header=[four_bytes(1:3);two_bytes(7:200)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ntr=get_no_of_traces(filename,nsamples)
% Use number of bytes to compute number of traces

    ll=dir(filename);
    nbytes=ll.bytes;
    ntr=0.25*(nbytes-3600)/(nsamples+60);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [headers]=read_trace_headers(fid,ntraces,nsamples)
% Read headers and traces
skip_trace = nsamples*4;
    headers = zeros(60,ntraces);
    for ii=1:ntraces
        hh4=fread(fid,60,'int32'); % read 240 bytes as four byte
        headers(:,ii) = hh4;
        %iloc=ftell(fid);
        fseek(fid,skip_trace,'cof');
%     fseek(fid,iloc,'bof');
%     hh2=fread(fid,120,'int16'); % read 240 bytes as two byte
%       if isempty(hh2)  % End of file reached
%          traces(:,ii:end)=[];
% 	 headers(:,ii:end)=[];
% 	 break
%       else
%          headers(true4four,ii)=hh4(indices(true4four));
%       end
%       headers(~true4four,ii)=hh2(indices(~true4four));
%       temp=fread(fid,param.no_samples,precision);
%       traces(:,ii)=temp(idx4times);
%       iloc=ftell(fid);
    if (floor(ii/20000) == ii/20000)
        fprintf('Scanning trace %d of %d (%d%% complete)\n',ii,ntraces,round((ii/ntraces)*100));
    end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [byte_value] = get_trace_header(n_traces,trace_header,index)   
        for ij=1:1:n_traces
            byte_value(ij) = trace_header(index,ij);
        end  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

