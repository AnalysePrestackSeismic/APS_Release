% Program Name
% Date: 25-April-2012
% Version: 0.0.1
% Authors: James Selvage
% What it does:
% Segy reader:
% Amplitude anomaly hunter:
% Structural closure finder:
% Constrained dix inversion:
% Intercept gradient inversion:
% Simultaneous prestack inversion:
% Constrained water column fitter:
% Phase estimation:
% MCMC deconvolution:

% Display all algorithms
[algo_name function_name] = process_to_run(1);

% Prompt depending on algorithm chosen

% Start editing here %
input_dir = {'/data/NOR/segy/knarr_2012_geotrace_prestm/post_prestm_residual_rms_vels'}; % add as many directories as needed
% paths to files instead?

% function call to scan all files in directories
files_in = directory_scan(input_dir);

% select files to input
%index_files = input('Enter numbers of segy files to scan (in bracket [], e.g. [1 3 5]): ');

switch function_name
    case 'vrms_vint_inversion'
        index_files = input('Enter number of Vrms segy file(s) to scan (in brackets [], e.g. [1 3 5]): ');
        % this will used for all files
        il_byte = 189;
        xl_byte = 193;
        scan_files = 1; % flag to scan all input files
        
        input_segy = segy_read_files(scan_files,index_files,files_in,il_byte,xl_byte);
        
        nfiles = length(index_files);
        
        % join together the scans
        for ii=1:1:nfiles
            trace_bytes_all_files{ii,1} = input_segy{ii}.trace_pointers;
        end
        
        trace_bytes_all_files = cell2mat(trace_bytes_all_files);
        
        aperture = input('Enter aperture: ');
        
        % enter il_inc and xl_inc
        il_inc = input_segy{1}.il_inc;
        xl_inc = input_segy{1}.xl_inc;
        
        % add filepath to array
        positions = segy_make_proc_pos_mult_files(il_inc,xl_inc,trace_bytes_all_files,aperture);
        
        % make array of file ids and paths
        for ii=1:1:nfiles
            filepaths_id = [input_segy{ii}.filepaths input_segy{ii}.file_id];
        end
        
        % Divide up input data
        max_cpus = 2; % since this is running as task it does not really matter
        % Control the memory usage
        n_blocks = 100;
        
        % Create matlab files for
        
    case 'sim_pre_stack_inversion'
        
end

% Scan all of the input files
% this will used for all files
%il_byte = 189;
%xl_byte = 193;
%extra_bytes_to_scan = [181 185];
% Do you want to scan all files that you have selected? [1 - Yes, 0 - No]
%scan_files = 1;
% processing with aperture?
%aperture = 1;
nfiles = length(index_files);

% scan the input files
%input_segy = segy_read_files(scan_files,index_files,files_in,il_byte,xl_byte,extra_bytes_to_scan);

% join together the scans
%for ii=1:1:nfiles
%    trace_bytes_all_files{ii} = input_segy{ii}.trace_pointers;
%end

%trace_bytes_all_files = cell2mat{trace_bytes_all_files};

% Maximum number of workers in cluster configuration
max_workers = 2; % since this is running as task it does not really matter
% Control the memory usage
n_blocks = 100;
output_dir = '/data/NOR/segy/knarr_2012_geotrace_prestm/matlab_out';

aperture = input('Enter aperture: ');

input_segy = segy_make_proc_pos(input_segy,aperture);
input_segy = njobs_worker(input_segy,max_workers,n_blocks);

[job tasks] = create_job_tasks(input_segy,function_name,output_dir,max_workers,n_blocks);

% Submit the job to be processed            
fprintf('Submitting the job to the cluster\n');
%submit(job);

% Get information about the job for the user
%get(job, 'Tasks')

% error checks on proc files etc.

% error checks on file types

% n_workers to submit to

% batch no wait so that we can have n_workers max and submit other jobs

