function [traces] = segy_read_traces(seismic,start_index,end_index,proc_flag,hor_flag)
% need to add range of loading

    seismic.fid = open_segy_file(seismic.filepaths);
    
    %traces = zeros(seismic.n_samples,(end_index-start_index)+1);
    
    if proc_flag == 0
        ii = 1;
        for k=start_index:1:end_index
            if isnan(seismic.trace_pointers(k,3))
                traces.data(:,ii) = zeros(seismic.n_samples,1);
            else
                fseek(seismic.fid,seismic.trace_pointers(k,3),'bof');
                temp = fread(seismic.fid,seismic.binary_header(8),'uint32=>uint32');
                traces.data(:,ii) = ibm2double(temp);
            end
            traces.pos(1,ii) = seismic.trace_pointers(k,1);
            traces.pos(2,ii) = seismic.trace_pointers(k,2);
            ii = ii+1;
        end
    
    elseif proc_flag == 1 && hor_flag == 0 % if has been calculated with an aperture
        ii = 1;
        for k=start_index:1:end_index
            if isnan(seismic.process(k,3))
                traces.data(:,ii) = zeros(seismic.n_samples,1);
            else
                fseek(seismic.fid,seismic.process(k,3),'bof');
                temp = fread(seismic.fid,seismic.n_samples,'uint32=>uint32');
                traces.data(:,ii) = ibm2double(temp);
            end
            traces.pos(1,ii) = seismic.process(k,1);
            traces.pos(2,ii) = seismic.process(k,2);
            ii = ii+1;
        end
    elseif proc_flag == 1 && hor_flag == 1    
        ii = 1;
        for k=start_index:1:end_index
            if isnan(seismic.process(k,4))
                traces.data(:,ii) = zeros(seismic.proc_samples,1);
            else
                fseek(seismic.fid,seismic.process(k,4),'bof');
                temp = fread(seismic.fid,seismic.proc_samples,'uint32=>uint32');
                traces.data(:,ii) = ibm2double(temp);
            end
            traces.pos(1,ii) = seismic.process(k,1);
            traces.pos(2,ii) = seismic.process(k,2);
            ii = ii+1;
        end        
    end
    
    traces.start_index = start_index;
    traces.end_index = end_index;
    
    close_segy_file(seismic.fid);
    
end

function data2=ibm2double(data)
% Convert IBM 32-bit floating-point format to IEEE 64-bit floating-point format
%      (based on the algorithm of function "ibm2num" written by Brian Farrelly)
% See also: ibm2single
%
% Written by: E. Rietsch: July 19, 2009 
% Last updated:
%
%         data2=ibm2double(data)
% INPUT 
% data    a matrix of unsigned 4-byte integers (uint32) representing data
%         in IBM floating-point format
% OUTPUT
% data2   corresponding matrix of doubles (IEEE format)

% get sign from first bit
% get exponent from first byte, last 7 bits
% remove bias from exponent 
% get mantissa from last 3 bytes

data2=(1-2*double(bitget(data,32))).*16.^ ...
  (double(bitshift(bitand(data,uint32(hex2dec('7f000000'))),-24))-64) .* ...
  (double(bitand(data,uint32(hex2dec('00ffffff'))))/2^24);
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