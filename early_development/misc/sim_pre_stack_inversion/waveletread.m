function [wavelet_data] = waveletread(input_wavelets_dir)
% Figures out the number of files in the current directory and reads the
% file names into a cell array. Then reads each file as ascii into a cell
% array and outputs. The format of the wavelet files must be in columns of
% TIME(ms), AMPLITUDE and without any header. The first two characters of the
% wavelet filename must be the angle of the stack it has come from (e.g.
% 05_wavelet_nears or 35_wavelet_fars).

start_point = pwd;
cd(input_wavelets_dir)

% Figure out the number of files in the current directory
[~,nfiles] = (system('ls -B | wc -l')); 
nfiles=str2double(nfiles);

% Read all filenames and convert from ascii to double
[~,fnames] = system('ls -B1');
numeric = double(fnames);

% Preallocate memory for some variables
count = 2;
fname_index = zeros(1,nfiles+1);
wavelet_data_cells{nfiles,1}=[];

% Loop to separate out each file name from one long character string
for i = 1:length(fnames)
    if numeric(1,i) == 10
        fname_index(1,count) = i;
        count = count+1;
    end
end

% Loop to read each file as ascii into a cell array
for i = 1:nfiles
    file_in = fnames(1,(fname_index(1,i)+1):(fname_index(1,i+1)-1));
    angle = str2double(file_in(1,1:2));
    if isnan(angle)
        error('Wavelet filenames must begin with the angle of the stack it has come from e.g ''05_wavelet'' or ''25_wavelet''!')
    end
    wavelet_data_cells{i,1} = [0,angle;dlmread(file_in)];
    wavelet_data(:,i+1) = wavelet_data_cells{i,1}(:,2);
end
wavelet_data(:,1) = wavelet_data_cells{1,1}(:,1);

wavelet_data = sortrows(wavelet_data',1)';

% Go back to the directory you were in to start with
cd(start_point);
end

