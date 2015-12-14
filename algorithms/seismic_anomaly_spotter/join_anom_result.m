function join_anom_result(job_meta_path,filepath,filename_string)
%% Defination: Joins the SAS result of all the blocks together, creates input for anomalous_body_connector

% Inputs:
%   job_meta_path: Path to job meta file
%   filepath: This is the folder name for reading the inputs and  writing the outputs into
%   filename_string: This is the search string with which the program will
%   scan the directory for input files

% Outpus:
%   returns no output
%   Writes on disk ...
%% Loading and processing inputs

job_meta = load(job_meta_path);                                 % Load job meta file

[files_in,nfiles] = directory_scan(filepath,filename_string);   % Scan firectory specicied in filepatn with search string fileame_string
files_in.names = sort_nat(files_in.names);                      % Natural sort file names
file_parts = regexp(files_in.names, '\_', 'split');

start_point = pwd;                                              % remember starting directory
filenames = '';                                                 % Intialize string for writing file names ??


% run using a little perl threading program
for i_file = 1:1:nfiles
    n_traces(i_file) = file_parts{i_file}(9); 
    n_traces(i_file) = regexprep(n_traces(i_file), 'traces', '');    
    ii_blocks{i_file} = str2double(file_parts{i_file}(2));   

    
end
samp = regexp(file_parts{end}(8),'\-','split');
samp = cell2mat(samp{1}(end));
%samp = cell2mat(samp{1});
ii_blocks = cell2mat(ii_blocks);
ii_blocks = sort(unique(ii_blocks));
n_blocks = length(ii_blocks);
n_parts = nfiles / n_blocks;

count = 1;
str_date = date;
str_date = regexprep(str_date, '-', '');

%% Part where the bits are glued together
file_out_bin = [filepath{1} 'sas_combine_sil' num2str(job_meta.pkey_min) 'sxl' num2str(job_meta.skey_min) 'nxl-' num2str(((job_meta.skey_max-job_meta.skey_min)/job_meta.skey_inc)+1) '-' 'samp' samp '_' str_date '.bin'];

fid_write = fopen(file_out_bin,'w');
for i_block = 1:1:n_blocks
    for i_part = 1:1:n_parts
        filepath_file = files_in.path{1};
        fid_read = fopen([filepath_file files_in.names{count}],'r');
        dataprob{i_part} = fread(fid_read,'float32');
        fclose(fid_read);
      
        dataprob{i_part} = reshape(dataprob{i_part},str2double(n_traces{count}),[]);
        count = count + 1;
    end
    dataprob = cell2mat(dataprob)';
    
   fwrite(fid_write,dataprob,'float32'); 
    
    clear dataprob
    
    % need to add unflatten
    fprintf('Completed block %d of %d\n',i_block,n_blocks);

end
fclose(fid_write);
%% Updating the job meta file 

job_meta.anom_result = file_out_bin;
save(job_meta_path,'-struct','job_meta','-v7.3');

end