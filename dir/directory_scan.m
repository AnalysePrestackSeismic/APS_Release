function [filter_files,nfiles] = directory_scan(input_dir,filename_string)
%% ------------------ Disclaimer  ------------------
% 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) makes no representation or warranty, express or implied, in 
% respect to the quality, accuracy or usefulness of this repository. The code
% is this repository is supplied with the explicit understanding and 
% agreement of recipient that any action taken or expenditure made by 
% recipient based on its examination, evaluation, interpretation or use is 
% at its own risk and responsibility.
% 
% No representation or warranty, express or implied, is or will be made in 
% relation to the accuracy or completeness of the information in this 
% repository and no responsibility or liability is or will be accepted by 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) in relation to it.
%% ------------------ License  ------------------ 
% GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
%% github
% https://github.com/AnalysePrestackSeismic/
%% ------------------ FUNCTION DEFINITION ---------------------------------
% compile_function: function to compile a matlab function as C++ code.
% Requires the MATLAB Compiler Toolbox. Algorithms should be stored in a
% subdirectory of the algorithms folder with the same name as the function
% to be compiled.
%   Arguments:
%       input_dir = cell array Input Directories
%       filename_string = the search string used a wildcard
%
%   Outputs:
%       filter_files = strucutre of file names and file paths of files 
%       found in the directory with the search string
%       nfiles = Number of filtered files found
%
%   Writes to Disk:
%       nothing

%%
n_dir = max(size(input_dir));   % Number of directories
start_point = pwd;              % Remember start directory and cd to input_dir directory

% -----------------SCAN EACH DIRECTORY----------------------------------
cumm_files = 1;
for ii=1:1:n_dir
    cd(input_dir{ii});          % change directory to iith lthe input directory
    
    % Figure out the number of files in the current directory
    [~,nfiles] = (system('ls -B * | wc -l')); 
    nfiles=str2double(nfiles);
    
    % Read all filenames and convert from ascii to double
    [~,fnames] = system('ls -B1 *');
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
     
    cumm_files = cumm_files + nfiles; % Increment cumm_files by number of files found in current directory  
    
end
files_in.nfiles = nfiles;
%-------------------------------------------------------------------------

%--------------FILTER THE FILES WITH THE SEARCH STRING-------------------

% fprintf('\nThe following files have been found:\n')
% for il = 1:files_in.nfiles
%     fprintf('File %d: %s in %s\n',il,files_in.names{il},files_in.path{il});
% end
%fprintf('\nApplying string filter, %s:\n',filename_string)

ii = 1;
for il = 1:files_in.nfiles
    if strfind(files_in.names{il},filename_string)          % if the file has the search string
%         fprintf('File %d: %s in %s\n',il,files_in.names{il},files_in.path{il});
        filter_files.names{ii} = files_in.names{il};        % Read file name in the structure
        filter_files.path{ii} = files_in.path{il};          % Read the file path in the structure
        ii = ii + 1;                                        % Increment the index of filtered file to be found
    end
end
nfiles = ii-1;                                              % Number of filtered files
%------------------------------------------------------------------------

cd(start_point);                                            % Go back to the directory you were in to start with

end