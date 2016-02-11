function [ ] = mat_to_segy(mode,hdr_mode,traces,headers,srate,ebcdic,outfile)
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
%MAT_TO_SEGY Writes matrices of traces and headers to a segy file
%   'mode':     keep/overwrite/append (default overwrite)
%   'hdr_mode': none/min/full (default min) 
%   traces:     nsamples x ntraces matrix of samples (floating point)
%   headers:    if hdr_mode=min then supply 3 x ntraces of inline, crossline, offset values
%               if hdr_mode=full then supply 60 x ntraces of header values
%               if hdr_mode=none then this parameter is ignored
%               all headers will be converted to integer if not already
%   srate:      sample rate in microseconds
%   ebcdic:     text string up to max size of 3200 characters
%   'outfile':  output segy file
%
%%=========================================================================
%   set up some variables
%
ntraces = size(traces,2);
nsamples = size(traces,1);
file_exists = 0; % flag for append mode
%%=========================================================================
%   set up header array if not supplied
%
if strcmp(hdr_mode,'none')
    headers = zeros(ntraces,60,'int32');
    headers(1,:) = 1:ntraces;
    headers(2,:) = 1:ntraces; % write trace num into 1st and 2nd header
    headers(8,:) = 65536; % writes 1 as 16-bit into first half of header
end

%%=========================================================================
%   check inputs
%
% if it looks like the headers are supplied with swapped rows and columns
% then transpose them
if or(size(headers,2)==3,size(headers,2)==60) && size(headers,1)==ntraces
    disp('Transposing header array');
    headers = headers';
end

if size(headers,2) ~= ntraces
    error('Number of headers and traces do not match');
end

if strcmp(hdr_mode,'min') && size(headers,1) ~= 3
        error('For hdr_mode=min you must supply three headers for each trace');
end

if strcmp(hdr_mode,'full') && size(headers,1) ~= 60
    error('For hdr_mode=min you must supply 60 headers for each trace');
end

if  ~isa(traces,'single')
    disp('Trace array is not single precision floating point and will be converted');
    traces = single(traces);
end

if ~isa(headers,'int32')
    disp('Header array is not 32 bit signed integer and will be converted');
    headers = int32(headers);
end

[fid, errmess] = fopen(outfile,'r'); % test if file exists by trying to read it

if fid ~= -1 % if file exists
    
    if strcmp(mode,'keep')
        if fid ~= -1
            error('Output file already exists');
        end
    else if strcmp(mode,'append')
            file_exists = 1;
            % do some more checks here that file is same srate etc as
            % incoming traces
            fseek(fid,3200,'bof');
            bin_header_in = fread(fid,400,'uint8');     % skip ebcdic and read the Binary Header
            
            % Re-interpret binary header as uint16 or uint32 as required
            two_bytes = bin_header_in(1:2:399)*256 + bin_header_in(2:2:400);
            four_bytes = ((bin_header_in(1:4:9)*256 + bin_header_in(2:4:10))*256+bin_header_in(3:4:11))*256+bin_header_in(4:4:12);
            bin_header_in = [four_bytes(1:3);two_bytes(7:200)];
            
            if bin_header_in(6) ~= srate
                error('Sample rate does not match existing file.');
            end
            
            if bin_header_in(8) ~= nsamples
                error('Number of samples does not match existing file.');
            end
            
        end
    end
    
    fclose(fid);
end


if size(ebcdic,2) > 3200 
    disp('ebcdic header string greater than 3200 characters and will be truncated');
end


%%=========================================================================
%   build header array if mode=min
% 
if strcmp(hdr_mode,'min')
    headers_out = zeros(60,ntraces,'int32');
    headers_out(1,:) = 1:ntraces; % trace number
    headers_out(2,:) = 1:ntraces; % trace number
    headers_out(8,:) = 65536; % writes 1 as 16-bit into first half of trace_id header
    headers_out(10,:) = headers(3,:); % offset
    headers_out(48,:) = headers(1,:); % inline
    headers_out(49,:) = headers(2,:); % crossline
    headers = headers_out;
    clear headers_out;
end

%%=========================================================================
%   pad or truncate ebcdic header to 3200 characters

ebcdic = sprintf('%-3200.3200s',ebcdic); 

%%=========================================================================
%   set up binary header
bin_header = zeros(200,1,'uint16'); % Define zero vector as size of trace header
bin_header(7) = 1;                  % Number of data traces per ensemble
bin_header(9) = srate;              % Sample rate (needs to be in microseconds)
bin_header(11) = nsamples;          % Number of samples per trace
bin_header(13) = 5;                 % SEGY data number format 5 = ieee float 754
bin_header(14) = 1;                 % Ensemble fold
bin_header(15) = 4;                 % Trace sorting code 4 = horizontally stacked , 2 = cdp ensemble
bin_header(28) = 1;                 % Measurement system
bin_header(151) = 256;              % SEGY revision number
bin_header(152) = 1;                % Fixed length flag  1= size always the same
       
%%=========================================================================
%   open and write binary and ebcdic headers
%
disp(['Opening ',outfile,' in ',mode,' mode']);
%
if strcmp(mode,'append') 
    fid = fopen(outfile,'a');
else
    fid = fopen(outfile,'w');
end

if ~strcmp(mode,'append') || ~file_exists % write binary and ebcdic if we're writing to a new file
    
    fwrite(fid,ebcdic,'char*1',0,'ieee-be'); % write ebcdic
    fwrite(fid,bin_header,'uint16',0,'ieee-be'); % write binary
    
end

%%=========================================================================
%   write traces and headers

traces = typecast(single(reshape(traces,1,(nsamples*ntraces))),'int32'); % typecast so matlab thinks traces are int32
segy_traces = [headers; reshape(traces,nsamples,ntraces)]; % then vertically concatenate with headers to make traces to write
    
fwrite(fid,segy_traces,'int32',0,'ieee-be'); %and write out as big endian ieee segy

disp(['Finished writing to ',outfile,' Closing file']);

fclose(fid);

end

