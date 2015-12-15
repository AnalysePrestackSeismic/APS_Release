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
[algo_name func_process] = process_to_run(1);

% Prompt depending on algorithm chosen

% Start editing here %
input_dir = {'/data/KEN/segy/2012_l10ab_geotrace/final_velocities_time/residual_rms_vel'}; % add as many directories as needed no slash at end
% paths to files instead?

% function call to scan all files in directories
files_in = directory_scan(input_dir);

% select files to input
index_files = input('Enter numbers of segy files to scan (in bracket [], e.g. [1 3 5]): ');
% this will used for all files
il_byte = 189;
xl_byte = 193;
extra_bytes_to_scan = [181 185];
% Do you want to scan all files that you have selected? [1 - Yes, 0 - No]
scan_files = 1;
% processing with aperture?
%aperture = 1;
nfiles = length(index_files);

% scan the input files
input_segy = segy_read_files(scan_files,index_files,files_in,il_byte,xl_byte,extra_bytes_to_scan);

% create processing positions file
%input_segy = segy_make_proc_pos(input_segy,aperture); % horizon_flag - should these be horizon bytes and position instead?

% % Maximum number of workers in cluster configuration
% max_workers = 20; 
% % Control the memory usage
% n_blocks = 1;
% output_dir = '/lustrecache/NOR/matlab_out';
% 
% %process = 'dix';
% %process = 'mcmc_decon';
% 
% % file_to_process = input('Enter index of file to process: ');
% aperture = 1;
% 
% %input('Enter aperture: ');
% % 
% input_segy = segy_make_proc_pos(input_segy,aperture);
% input_segy = njobs_worker(input_segy,max_workers,n_blocks);
% 
% 
% tol = 1e-6;
% iter = 25;
% 
% output_location = output_dir;
% 
%         for ii = 1:1:max_workers*n_blocks
% %             %tasks(ii) = createTask(job,fnc_handle,0,{ ...          
% %             %    input_files{file_to_process}.filepaths, ...
% %             %    input_files{file_to_process}.aperture, ...
% %                 input_files{file_to_process}.n_samples, ...
% %                 input_files{file_to_process}.process(input_files{1}.index_worker(ii,4):input_files{1}.index_worker(ii,5),:), ...
% %                 tol, ...
% %                 iter, ...
% %                 output_dir});
%             traces_process = input_segy{1}.process(input_segy{1}.index_worker(ii,4):input_segy{1}.index_worker(ii,5),:);
%             filepaths = input_segy{1}.filepaths;
%             aperture = input_segy{1}.aperture;
%             n_samples = input_segy{1}.n_samples; 
%   % filepaths,aperture,n_samples,traces_process,tol,iter,output_location          
%             fileid = num2str(ii);
%             job = strcat('job',fileid);
%             file_mat = strcat('job',fileid,'.mat');
%             save(file_mat,'filepaths','job','aperture','n_samples','traces_process','tol','iter','output_location');
%         end
% 
% 
% % 
% % [job tasks] = create_job_tasks(input_segy,file_to_process,func_process,output_dir,max_workers,n_blocks);
% % 
% % % Submit the job to be processed            
% % fprintf('Submitting the job to the cluster\n');
% %submit(job);
% 
% % Get information about the job for the user
% %get(job, 'Tasks')
% 
% % error checks on proc files etc.
% 
% % error checks on file types
% 
% % n_workers to submit to
% 
% % batch no wait so that we can have n_workers max and submit other jobs

