function [files_in] = directory_scan(input_dir)

% Number of directories
n_dir = max(size(input_dir));

% Remember start directory and cd to the input_logs directory
start_point = pwd;

cumm_files = 1;
% Scan each directory
for ii=1:1:n_dir
    cd(input_dir{ii});
    
    % Figure out the number of files in the current directory
    [~,nfiles] = (system('ls -B *gy | wc -l')); 
    nfiles=str2double(nfiles);
    
    % Read all filenames and convert from ascii to double
    [~,fnames] = system('ls -B1 *gy');
    numeric = double(fnames);
    
    % Preallocate memory for some variables
    count = 2;
    fname_index = zeros(1,nfiles+1);
    
    % Loop to separate out each file name from one long character string
    for ij= 1:length(fnames)
        if numeric(1,ij) == 10
            fname_index(1,count) = ij;
            count = count+1;
        end
    end
    
    % Loop to read each file as ascii into a cell array
    for ik=1:nfiles
        files_in.names{cumm_files+ik-1} = fnames(1,(fname_index(1,ik)+1):(fname_index(1,ik+1)-1));
        files_in.path{cumm_files+ik-1} = input_dir{ii};
    end
     
    cumm_files = cumm_files + nfiles;   
    
end

files_in.nfiles = nfiles;

fprintf('\nThe following files have been found:\n')
for il = 1:files_in.nfiles
    fprintf('File %d: %s in %s\n',il,files_in.names{il},files_in.path{il});
end

% Go back to the directory you were in to start with
cd(start_point);

end