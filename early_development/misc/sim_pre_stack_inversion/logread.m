function [log_data] = logread(input_logs_dir)
% Figures out the number of files in the current directory and reads the
% file names into a cell array. Then reads each file as ascii into a cell
% array and outputs. The format of the log files must be in columns of TIME(ms),
% AI, SI, RHO and without and header. The first 8 characters in the
% filename must be the INLINE location of the well. The second 8 characters
% in the filename must be the CROSSLINE location of the well. For example:
% 0000123400005678_well_log_file.txt would mean the well is at inline 1234
% and crossline 5678.

% Remember start directory and cd to the input_logs directory
start_point = pwd;
cd(input_logs_dir);

% Figure out the number of files in the current directory
[~,nfiles] = (system('ls -B | wc -l')); 
nfiles=str2double(nfiles);

% Read all filenames and convert from ascii to double
[~,fnames] = system('ls -B1');
numeric = double(fnames);

% Preallocate memory for some variables
count = 2;
fname_index = zeros(1,nfiles+1);
log_data{nfiles,1}=[];

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
    il = str2double(file_in(1,1:8));
    xl = str2double(file_in(1,9:16));
    if (isnan(il)||isnan(xl))
        error('Well log filenames must begin with the inline and crossline location of the well e.g ''0000123400005678_well''=il1234 xl5678!')
    end
    log_data{i,1} = [il,xl,0,0;dlmread(file_in)];
    log_data{i,1}(1,3) = min(log_data{i,1}(2:end,1));
    log_data{i,1}(1,4) = max(log_data{i,1}(2:end,1));
end

% Go back to the directory you were in to start with
cd(start_point);
end

