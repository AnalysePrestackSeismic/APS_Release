function [] = optslice(file_in,bytespersample,mem)
% Optimises opendtect binary files for time slice access.
%
% USAGE:   [] = optslice(file_in,bytespersample,mem)
%
% INPUTS:   file_in = binary file from opendtect with sampling info and inline and crossline locations in the file.
%           bytespersample = number of bytes per sample (4 if input data is float32 format).
%           mem = is the memory block in Gb for use by the code. Actual memory used will be higher than this.
%
% OUTPUTS:  sliceopt_file_in = slice optimised binary file (without sampling info or inline and crossline loacions).
           
% Open the input binary file.
fid1 = fopen(file_in);

% Create the output slice access optimised file.
fid2 = fopen(sprintf('sliceopt_%s',file_in),'a');

% Work out the size of the binary file in bytes.
fseek(fid1, 0, 'eof');
filesize = ftell(fid1);

% Read the number of time samples from binary file.
fseek(fid1,2*bytespersample,'bof');
nt = fread(fid1,1,'int');

% Change type from int to double to avoid divide errors.
nt = double(nt);

% Calculate the number of traces in the file.
ntrace = (filesize-3*bytespersample)/((nt+2)*bytespersample);

% Work out the block size in time samples to read such that the size of the block is equal to 'mem'.
block = mem*(1024^3)/(ntrace*bytespersample);

% Ensure the block size is smaller than the number of time samples.
if block >= nt
    block = nt;
end

% Calculate the number of loops needed to read the whole file.
loops = floor(nt/block);
block = floor(block);

% Pre-allocate memory for temporary variable to hold a block of data.
tmpslice = zeros(ntrace,block);

% Skip past the next two samples (first inline and crossline position) to the first data byte.
fseek(fid1,2*bytespersample,'cof');

% Loops to read and write the data.
for j = 1:loops
    % Print read progress to screen.
    fprintf('reading slices %d to %d of %d\n',1+((j-1)*block),j*block,nt)
    % Loop over all traces.
    for k = 1:ntrace
        % Read block of time samples on a trace.
        tmpil = fread(fid1,block,'*float32');
        % Jump to same block of time samples on the next trace.
        fseek(fid1,(2+nt-block)*bytespersample,'cof');
        % Build matrix of slices read.
        tmpslice(k,:) = tmpil';
    end
    % Print write progress to screen.
    fprintf('writing slices %d to %d of %d\n',1+((j-1)*block),j*block,nt)
    % Loop over all time samples in the block.
    for k = 1:block
        % Write time samples in block to file.
        fwrite(fid2,tmpslice(:,k),'float32');
    end
    % Jump to next block of time samples to read.
    fseek(fid1,(5+(j*block))*bytespersample,'bof');
end

% Loops to read the remaining part of the data.
if block*loops ~= nt
    leftovers = nt - (block*loops);
    tmpslice = zeros(ntrace,leftovers);
    fprintf('reading slices %d to %d of %d\n',j*block+1,nt,nt)
    for k = 1:ntrace
        tmpil = fread(fid1,leftovers,'*float32');
        fseek(fid1,(2+nt-leftovers)*bytespersample,'cof');
        tmpslice(k,:) = tmpil';
    end
    fprintf('writing slices %d to %d of %d\n',j*block+1,nt,nt)
    for k = 1:leftovers
        fwrite(fid2,tmpslice(:,k),'float32');
    end
end

% Close all files.
fclose all;

end