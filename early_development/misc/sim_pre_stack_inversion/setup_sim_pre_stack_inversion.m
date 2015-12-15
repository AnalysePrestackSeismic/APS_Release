% Program Name
% Date: 25-April-2012
% Version: 0.0.1
% Authors: James Selvage
% What it does:

% Start editing here
input_dir = {'/lustrecache/TZA/segy/2012_kussini_pgs_full_sequence/angle_stacks'}; % add as many directories as needed no slash at end
output_location = '/lustrecache/TZA/matlab_out';
% paths to files instead?

input_logs_dir = '/lustrecache/TZA/dtect/kussini_2012_full_seq_pgs/Misc/matlab_inversion/input_logs';
log_data = logread(input_logs_dir);

input_wavelets_dir = '/lustrecache/TZA/dtect/kussini_2012_full_seq_pgs/Misc/matlab_inversion/input_wavelets';
wavelet_data = waveletread(input_wavelets_dir);

% Maximum number of workers in cluster configuration
max_workers = 20; 

% Control the memory usage
n_blocks = 1;

% Inversion settings
tol = 1e-6;
iter = 25;
k = 1.37602;
m = 0.18682;
kc = -4.02995;
mc = -0.796873;
gam = 0.4964;

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

% Use a horizon
input_hor  = input('Do you want to use a horizon? [1 - Yes, 0 - No]: ');

if (input_hor == 1)
    hor_path = input('Enter path to horizon file: ');
    seismic.horpath = hor_path;
     
    raw = dlmread(seismic.horpath, '\t');
    hor.data = sortrows(raw);
    clear raw;
    
    start_time = input('Enter start time (ms) (+ve for below horizon, -ve for above horizon): ');
    time_window = input('Time window to process (ms): ');
    
    for ii=1:1:nfiles
        hor.data = unique(hor.data,'rows','first');
        [hor.Lia, hor.Locb] = ...
            ismember(input_segy{ii}.trace_pointers(:,1:2),hor.data(:,1:2),'rows');
       
        window_length = (start_time/(input_segy{ii}.s_rate/1000));
        input_segy{ii}.trace_pointers(hor.Lia,4) = input_segy{ii}.trace_pointers(hor.Lia,3) + round(((hor.data(nonzeros(hor.Locb),3)/(input_segy{ii}.s_rate/1000))+1))*4 + window_length*4;

        input_segy{ii}.trace_pointers(~hor.Lia,4) = NaN;
    
        input_segy{ii}.proc_samples = 1+(time_window/(input_segy{ii}.s_rate/1000));
    end
end

% process without aperture
input_segy = segy_make_proc_pos(input_segy,0);
% create jobs
input_segy = njobs_worker(input_segy,max_workers,n_blocks);

% Work out number of angle stacks
ii = 1;
for i_file = 1:1:nfiles
    if input_segy{i_file}.type == 2
        angles(ii) = input_segy{i_file}.angle*pi()/180;
        angle_index(ii) = i_file;
        ii = ii + 1;        
    end
end

% Scale wavelets
% get inline crossline of log data
for i_well = 1:1:size(log_data,1)
   il_xl_pos(i_well,1:2) = [log_data{i_well,1}(1,1) log_data{i_well,1}(1,2)];
   % need this to row per log and 2 cols   
end

for i_angle=1:1:length(angle_index) % for each angle stack
    % make sure only angle stacks
    [Lia, Locb] = ismember(il_xl_pos(:,1:2),input_segy{angle_index(i_angle)}.trace_pointers(:,1:2), 'rows');
    %log_seis{ij}.traces(:,3) = NaN(length(input_segy{ij}.trace_pointers(:,1)),1);
    log_seis{i_angle}.trace_pointers(Lia,1) = input_segy{angle_index(i_angle)}.trace_pointers(nonzeros(Locb),1); 
    log_seis{i_angle}.trace_pointers(Lia,2) = input_segy{angle_index(i_angle)}.trace_pointers(nonzeros(Locb),2);
    log_seis{i_angle}.trace_pointers(Lia,3) = input_segy{angle_index(i_angle)}.trace_pointers(nonzeros(Locb),3);  
    log_seis{i_angle}.filepath = input_segy{angle_index(i_angle)}.filepaths;
    % get seismic data
    % log_seis{i_angle} = segy_read_traces(input_segy{angle_index(i_angle)}.filepaths,input_segy{angle_index(i_angle)}.n_samples,log_seis{i_angle}.trace_pointers); 
end

clear Lia Locb

% calculate coefficients for inversion and wavelet scaling
[c1 c2 c3] = calculate_coefficients(angles,0,gam,k,m);

for i_well = 1:1:size(log_data,1)
% scale wavelets
    
    ref_imp = zeros((length(log_data{i_well,1}(:,1)))-1,3);
    for l = 2:4
        ref_imp(:,l-1) = [diff(log_data{i_well,1}(2:end,l));0]./conv2(log_data{i_well,1}(2:end,l),[1,1],'same');   
    end
    ref_imp = ref_imp(1:end-1,:);

    for i_angle=1:1:length(angle_index)
        % Calculate angle dependent reflectivity and accompanying synthetic
        % seismograms
        ref_data{i_well,1}(:,i_angle) = [angles(i_angle);(c1(i_angle,1)*ref_imp(:,1))+(c2(i_angle,1)*ref_imp(:,2))+(c3(i_angle,1)*ref_imp(:,3))];
        syn_data{i_well,1}(:,i_angle) = [angles(i_angle);conv2(ref_data{i_well,1}(2:end,i_angle),wavelet_data(2:end,i_angle+1),'same')];
        trace = segy_read_traces(input_segy{angle_index(i_angle)}.filepaths,input_segy{angle_index(i_angle)}.n_samples,log_seis{i_angle}.trace_pointers(i_well,:));
        trace_data{i_well,1}(:,i_angle) = [angles(i_angle);trace.data];
        scalar(i_well,i_angle) = sqrt(sum(trace_data{i_well,1}(2:end,i_angle).^2))/(sqrt(sum(ref_data{i_well,1}(2:end,i_angle).^2))*sqrt(sum(wavelet_data(2:end,i_angle+1).^2)));
    end
    
end

avg_scalar = mean(scalar,1);
scaled_wavelet_data(:,1) = wavelet_data(:,1);
scaled_wavelet_data(1,:) = wavelet_data(1,:);

for i_angle=1:1:length(angle_index)
    scaled_wavelet_data(2:end,i_angle+1) = wavelet_data(2:end,i_angle+1)*avg_scalar(1,i_angle);
end




% Get background model


% Make job files
for ii = 1:1:max_workers*n_blocks
    traces_process = input_segy{1}.process(input_segy{1}.index_worker(ii,4):input_segy{1}.index_worker(ii,5),:);
    filepaths = input_segy{1}.filepaths;
    aperture = input_segy{1}.aperture;
    n_samples = input_segy{1}.n_samples; 
    fileid = num2str(ii);
    job = strcat('job',fileid);
    file_mat = strcat('job',fileid,'.mat');
    save(file_mat,'filepaths','job','aperture','n_samples','traces_process','tol','iter','output_location');

    
end

    
%     
%     % number of angle stacks
%    for 
%    
%    % number of traces into cube
% 
% 
%        
%      trace_byte =    
%        
%        
%    input_segy{ii}.min_inline
%    input_segy{ii}.seismic.min_xline
%    
%    extract_trace = 
%     
% end



% for ii = 1:1:max_workers*n_blocks
% %             %tasks(ii) = createTask(job,fnc_handle,0,{ ...          
% %             %    input_files{file_to_process}.filepaths, ...
% %             %    input_files{file_to_process}.aperture, ...
% %                 input_files{file_to_process}.n_samples, ...
% %                 input_files{file_to_process}.process(input_files{1}.index_worker(ii,4):input_files{1}.index_worker(ii,5),:), ...
% %                 tol, ...
% %                 iter, ...
% %                 output_dir});
%     for ij = 1:1:length(angles)
%             traces_process{ij} = input_segy{ij}.process(input_segy{ij}.index_worker(ii,4):input_segy{ij}.index_worker(ii,5),:);
%             filepaths{ij} = input_segy{ij}.filepaths;
%             
%             n_samples = input_segy{1}.n_samples; 
%   % filepaths,aperture,n_samples,traces_process,tol,iter,output_location          
%             fileid = num2str(ii);
%             job = strcat('job',fileid);
%             file_mat = strcat('job',fileid,'.mat');
%             save(file_mat,'filepaths','job','aperture','n_samples','traces_process','tol','iter','output_location'); angles
%     end
%         enmd


% 
% [job tasks] = create_job_tasks(input_segy,file_to_process,func_process,output_dir,max_workers,n_blocks);
% 
% % Submit the job to be processed            
% fprintf('Submitting the job to the cluster\n');
%submit(job);

% Get information about the job for the user
%get(job, 'Tasks')

% error checks on proc files etc.

% error checks on file types

% n_workers to submit to

% batch no wait so that we can have n_workers max and submit other jobs

