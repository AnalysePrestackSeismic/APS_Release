function [] = geomstrip(file_in,bytespersample)
% Strips geometry info from an opendtect binary file containing sampling info and inline and crossline locations.
%
% USAGE:   [] = geomstrip(file_in,bytespersample)
%
% INPUTS:   file_in = opendtect binary file.
%           bytespersample = number of bytes per sample (4 if input data is float32 format).
%
% OUTPUTS:  geometry_file_in = geometry file containing sampling info and inline and crossline loacions for file_in. This can be used as input to the function optil.

% Open input file.
fid1 = fopen(file_in);

% Create output geometry file.
fid2 = fopen(sprintf('geometry_%s',file_in),'a');

% Work out the size of the input file.
fseek(fid1, 0, 'eof');
filesize = ftell(fid1);
frewind(fid1);

% Write start time to geometry file.
fwrite(fid2,fread(fid1,1,'*float32'),'float32');

% Write sample rate to geometry file.
fwrite(fid2,fread(fid1,1,'*float32'),'float32');

% Read the number of time samples from input file.
nt = fread(fid1,1,'*int');

% Write number of time samples to geometry file.
fwrite(fid2,nt,'int');

% Calculate the number of traces in the input file.
ntrace = (filesize-3*bytespersample)/((nt+2)*bytespersample);

% Initialise variable.
geometry = zeros(2*ntrace,1);

% Update user on progress.
fprintf('Reading geometry from file...\n');

% Read the geometry info from input file.
for k = 1:ntrace
    geometry(((2*k)-1):(2*k),1) = fread(fid1,[2,1],'2*int',nt*bytespersample);
end

% Update user on progress.
fprintf('Writing geometry to file...\n');

% Write the geometry info to output geometry file.
fwrite(fid2,geometry,'int');

% Close all files.
fclose all;

end
    

