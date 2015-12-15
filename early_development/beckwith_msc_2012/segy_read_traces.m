function [traces] = segy_read_traces(filepath,n_samples,traces_process)
% need to add range of loading

    fid = open_segy_file(filepath);
    
    n_traces = size(traces_process,1);
    traces.data = zeros(n_samples,n_traces);
    
    for ik=1:1:n_traces
        if isnan(traces_process(ik,3))
            traces.data(:,ik) = zeros(n_samples,1);
        else
            fseek(fid,traces_process(ik,3),'bof');
            temp = fread(fid,n_samples,'uint32=>uint32');
            traces.data(:,ik) = ibm2double(temp);
        end
            traces.pos(1,ik) = traces_process(ik,1);
            traces.pos(2,ik) = traces_process(ik,2);
    end
    
    close_segy_file(fid);    
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