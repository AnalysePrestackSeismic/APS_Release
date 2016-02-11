function[]=run_batch_2d_algo(input_filepath,varargin)
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
%% Function to run trim on a batch of 2D prestack gathers

% if you want scan and execution of an algo after that
% input_filepath = directory of input files as a cell array with input
% string. Dont end the path with '/'
% varargin : all scanning parameters as string , algorithm name, slurm_part, n_cores, then all algorithm
% parameters as strings


%if you dont want to scan, just want to run an algo in loop
% input_filepath = path of batch run info mat file, please give this as a
% string not as cell array
% varargin :  algorithm name, slurm_part, n_cores, then all algorithm
% parameters as strings

% prints : batch file run info file name
% saves to disk : bath file run info file apartfrom files associated with
% scanning and the specific algorithm
%%
f_srch=regexp(input_filepath,'_run_info','ONCE');           %Search if input file path has _run_info in it, if yes dont run scan
rg=1;   %initialize a counter for accessing varagin entries

%% ------------run the scan jobs and construct batch run info file ---------------------
if isempty(f_srch{1})
    %  Sort the variables form varargin
    filename_string=varargin{rg};rg=rg+1;
    output_filepath=varargin{rg};rg=rg+1;
    il_byte=varargin{rg};rg=rg+1;
    xl_byte=varargin{rg};rg=rg+1;
    offset_byte=varargin{rg};rg=rg+1;
    parallel=varargin{rg};rg=rg+1;
    anggath=varargin{rg};rg=rg+1;    
    algorithm_name=varargin{rg};rg=rg+1;
    
    
    % ----------------scan for files----------------
    [files_in,nfiles] = directory_scan(input_filepath,filename_string);     % files_in is a structure, .names cell array, .path cell array
    % directory_scan filters out the file names with the user supplied string and also returns the number of such files: nfiles that ...
    % ... can be used for index pre alocation in future steps
    %   files_in.names = sort_nat(files_in.names);                          % Natural Sort file names
    %start_point = pwd;                                                     % Remember starting directory
    
    % -------------make folders structure for input to scanning---------------------------------------------
   
    if exist(output_filepath,'dir') == 0
        mkdir(output_filepath)                                              % Make output directory (umbrella directory)
        system(['chmod 777 ',output_filepath]);                             % Give permissions
    end    
    links_dir=strcat(output_filepath,'links/');    
    if exist(links_dir,'dir') == 0
        mkdir(links_dir)                                                    % Make sub direcoty directory for storing individual link directories for all lines
        system(['chmod 777 ',links_dir]);                                   % Give permissions
    end    
    results_dir=strcat(output_filepath,'result/');                      
    if exist(results_dir,'dir') == 0                                        
        mkdir(results_dir)                                                  % Make output directory for storing results from algorithm(umbrella directory)
        system(['chmod 777 ',results_dir]);                                 % Give permissions
    end    
    % --------------------make links and output directories----------------------------------------------------
    % Loop through all the input files
    for f=1:nfiles
        file_name_temp=strcat(files_in.path{f},'/',files_in.names{f});      % Construct input file name
        id=strfind(files_in.names{f},'.');
        links_input_dir{f}=strcat(links_dir,files_in.names{f}(1:(id-1)),'/');% Construct soft link name 
        if exist(links_input_dir{f},'dir') == 0
            mkdir(links_input_dir{f})                                       % Make individual directory for storing soft links with in the link directory
            system(['chmod 777 ',links_input_dir{f}]);
        end
        link_name_temp=strcat(links_input_dir{f},'gather_',algorithm_name,'_input_link_block1.segy'); % Create string name for the soft link
        system(['ln -s ',file_name_temp,' ',link_name_temp]);               % Make the soft link
        output_dir{f}=strcat(results_dir,files_in.names{f}(1:(id-1)),'/');  % Make individual results directory name for storing the results for the lines
        
        if exist(output_dir{f},'dir') == 0
            mkdir(output_dir{f})                                            % Create results directory if it doesnt exist
            system(['chmod 777 ',output_dir{f}]);                           % Give permissions
        end
    end    
    clear f file_name_temp id  link_name_temp ;                             % Clear variables
    % ------------------- run segy_make job----------------------
    %Loop though all the input files
    for g=1:nfiles
        job_meta_files{g} = segy_make_job(links_input_dir{g},'block',il_byte,xl_byte,offset_byte,parallel,anggath,output_dir{g}); % scan the input gather with segy make job
        fprintf('Completed segy make job for %d out of %d lines \n',g,nfiles);                                                    % Print which file has been scanned
    end
    batch_run_file_info = struct;                                           % Define a structure for storing batch run information
    batch_run_info_path=strcat(output_filepath,'batch_run_info.mat');       % path of output batch run info file
    
    %tore various information about batch run in  batch run info file
    batch_run_file_info.files_in=files_in;
    batch_run_file_info.nfiles=nfiles;
    batch_run_file_info.links_input_dir=links_input_dir;
    batch_run_file_info.job_meta_files=job_meta_files;
    batch_run_file_info.il_byte=il_byte;
    batch_run_file_info.xl_byte=xl_byte;
    batch_run_file_info.offset_byte=offset_byte;
    batch_run_file_info.parallel=parallel;
    batch_run_file_info.anggath=anggath;
    
else
    algorithm_name=varargin{rg};rg=rg+1;
    batch_run_info_path=input_filepath{1};
    batch_run_file_info=load(batch_run_info_path);
end
%% ----------------------run algo--------------------------------------------

% if the job meta files exists

% slurm_part=varargin{rg};rg=rg+1;
% n_cores=varargin{rg};rg=rg+1;

% construct the tail string with the parameters if there are parametersleft
if rg<=length(varargin)
    param_script=[varargin{rg}];rg=rg+1;
    for kk=rg:length(varargin)
        param_script=[param_script,' ',varargin{kk},];
    end
else
    param_script=[];
end

if  (isempty(strfind(algorithm_name,'trim'))&&isempty(strfind(algorithm_name,'int_grad_inv_proj'))&&isempty(strfind(algorithm_name,'wavelet_estimation'))&&isempty(strfind(algorithm_name,'ava_analysis')))
    %for programs that dont run on cluster
    for g=1:batch_run_file_info.nfiles
        % for wavelet averaging
        if ~isempty((strfind(algorithm_name,'wavelet_avg')))
            text_script = ['/apps/gsc/matlab-library/development/maps/algorithms/wavelet_estimation/run_',algorithm_name,'.sh /apps/matlab/v2013a/',' ',batch_run_file_info.job_meta_files{g},' ',param_script];
        end
        % for smoothening meta data
        if ~isempty((strfind(algorithm_name,'meta_data_2d_smo_xy')))
            text_script = ['/apps/gsc/matlab-library/development/maps/extra/run_',algorithm_name,'.sh /apps/matlab/v2013a/',' ',batch_run_file_info.job_meta_files{g},' ',param_script];
        end
        system(text_script);
        fprintf('submitted line : %g \n',g)
    end
else
    % construct the run scripts and submit to the clutser
    for g=1:batch_run_file_info.nfiles
        %     text_script = ['/apps/gsc/matlab-library/development/maps/node/run_node_slurm_submit_neo.sh /apps/matlab/v2013a/',' ',algorithm_name,' ',batch_run_file_info.job_meta_files{g},' ',slurm_part,' ',n_cores,' ',param_script];
        text_script = ['/apps/gsc/matlab-library/development/maps/node/run_node_slurm_submit_neo.sh /apps/matlab/v2013a/',' ',algorithm_name,' ',batch_run_file_info.job_meta_files{g},' ',param_script];
        batch_run_file_info.last_algo_run_script{g}=text_script;    % Store the script that will be submitted to system.
        [~,text_return]=system(text_script);
        text_return2=deblank(regexp(text_return,'\n','split'));
        batch_run_file_info.algo_run_qc{g}=text_return2{1,9};
        pause(10);
        fprintf('submitted line : %g \n',g)
    end
end
%%
save(batch_run_info_path,'-struct','batch_run_file_info','-v7.3');             % Saves Seismic structure to mat file
disp(batch_run_info_path);                                                     % Display batch run infor file path


end