function seis = segy_read_files(file_path)

files_in = data_read(file_path);

fprintf('\nThe following files are in directory: %s\n',files_in.path)
for i = 1:files_in.nfiles
    fprintf('File %d: %s\n',i,files_in.names{i});
end

index_files = input('Enter number of segy file to scan (in bracket [], e.g. [1 3 5]): ');

il_byte = input('Enter inline byte: ');
xl_byte = input('Enter crossline byte: ');

scan_files = input('Do you want to scan all files that you have selected? [1 - Yes, 0 - No]: ');

if (scan_files == 1)
    for i = 1:length(index_files)
       seis{i} = segy_make_index(cell2mat(strcat(files_in.path,'/',files_in.names(index_files(i)))),il_byte,xl_byte);
    end
elseif (scan_files == 0)
    seis{1} = segy_make_index(cell2mat(strcat(files_in.path,'/',files_in.names(index_files(1)))),il_byte,xl_byte);   
    
    for i = 2:length(index_files)
        seis{i} = seis{1};
        seis{i}.filepaths = strcat(files_in.path,'/',cell2mat(files_in.names(index_files(i))));
    end
end

for i = 1:length(index_files)
    fprintf('\nFor file %s enter: \n',seis{i}.filepaths)
    seis{i}.type = input('Enter file type [1 - Full stack, 2 - Angle stack, 3 - Velocity, 4 - Attribute]: ');
        if  (seis{i}.type == 2)
            seis{i}.angle = input('Angle in degrees: ');
        end
        if  (seis{i}.type == 4)
            seis{i}.attribtype = input('Enter attribute type: ');
        end
end

end

function [files_in] = data_read(input_dir)

% Remember start directory and cd to the input_logs directory
start_point = pwd;
cd(input_dir);

% Figure out the number of files in the current directory
[~,nfiles] = (system('ls -B | wc -l')); 
nfiles=str2double(nfiles);

% Read all filenames and convert from ascii to double
[~,fnames] = system('ls -B1');
numeric = double(fnames);

% Preallocate memory for some variables
count = 2;
fname_index = zeros(1,nfiles+1);

% Loop to separate out each file name from one long character string
for i = 1:length(fnames)
    if numeric(1,i) == 10
        fname_index(1,count) = i;
        count = count+1;
    end
end

% Loop to read each file as ascii into a cell array
for i = 1:nfiles
    files_in.names{i} = fnames(1,(fname_index(1,i)+1):(fname_index(1,i+1)-1));
end

files_in.nfiles = nfiles;
files_in.path = input_dir;

% Go back to the directory you were in to start with
cd(start_point);
end