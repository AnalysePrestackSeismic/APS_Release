
%% Need to CD to right directory then run the following to end up with a cell array of filenames:


cumm_files = 1;
[~,nfiles] = (system('ls -B *gy | wc -l'));
nfiles=str2double(nfiles);

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
    names{cumm_files+ik-1} = fnames(1,(fname_index(1,ik)+1):(fname_index(1,ik+1)-1));
    % files_in.path{cumm_files+ik-1} = input_dir{ii};
end