function [job tasks indexes] = create_job_tasks(input_files,function_name,...
    output_dir,max_workers,n_blocks)

fnc_handle = str2func(function_name);

library_path = '/apps/gsc/matlab-mcode-beta/beta-library';

common_functions = '/apps/gsc/matlab-mcode-beta/beta-library/segy_tools/segy_read_traces.m';

file_to_process = 1;

tol = 1e-6;
iter = 25;

output_location = output_dir;

switch function_name
    case 'vrms_vint_inversion'
%         function_path = strcat(library_path,'/',function_name,'/',function_name,'.m');
%         
%         path_depend = {common_functions function_path};
%         
%         tol = input('Enter tolerance: ');
%         iter = input('Enter iterations: ');
% 
%         sched = findResource('scheduler','configuration', defaultParallelConfig);
% 
%         job = createJob(sched,'FileDependencies', ...
%                         path_depend);

tol = 1e-6;
iter = 25;

output_location = output_dir;

        for ii = 1:1:max_workers*n_blocks
%             %tasks(ii) = createTask(job,fnc_handle,0,{ ...          
%             %    input_files{file_to_process}.filepaths, ...
%             %    input_files{file_to_process}.aperture, ...
%                 input_files{file_to_process}.n_samples, ...
%                 input_files{file_to_process}.process(input_files{1}.index_worker(ii,4):input_files{1}.index_worker(ii,5),:), ...
%                 tol, ...
%                 iter, ...
%                 output_dir});
            traces_process = input_segy{1}.process(input_segy{1}.index_worker(ii,4):input_segy{1}.index_worker(ii,5),:);
            filepaths = input_segy{1}.filepaths;
            aperture = input_segy{1}.aperture;
            n_samples = input_segy{1}.n_samples; 
  % filepaths,aperture,n_samples,traces_process,tol,iter,output_location          
            fileid = num2str(ii);
            file_mat = strcat('job',fileid,'.mat');
            save(file_mat,'filepaths','aperture','n_samples','traces_process','tol','iter','output_location');
        end
    case 'sim_pre_stack_inversion'
        
        sched = findResource('scheduler','configuration', defaultParallelConfig);

        job = createJob(sched,'FileDependencies', ...
                        path_depend);
                
        for ii = 1:1:max_workers*n_blocks
            tasks(ii) = createTask(job,fnc_handle,0,{ ...          
                input_files{file_to_process}.filepaths, ...
                input_files{file_to_process}.aperture, ...
                input_files{file_to_process}.n_samples, ...
                input_files{file_to_process}.process(input_files{1}.index_worker(ii,4):input_files{1}.index_worker(ii,5),:), ...
                tol, ...
                iter, ...
                output_dir});
        end

end


% save('job1', )

% switch function_name
%     case sim_pre_stack_inversion
%         % number of files read in
%         n_files = length(input_files);
%         
%         % find angles from input files
%         ii = 1;
%         for i_file = 1:1:n_files
%                       
%            if input_files{i_file}.type == 2
%                angles(ii) = input_files{i_file}.angle;
%                ii = ii + 1;
%                nt(ii) = input_files{i_file}.n_samples;
%            end
%         end
%         
%         % ask user for background model
%         
%         % read / estimate wavelets
%         
%         % scale wavelets
%         
%         % number of angle stacks
%         n_angles = length(angles);
%         
%         no_rho = input('Do you want to invert for rho [1 - Yes, 0 - No]');
%         
%         [c1 c2 c3] = calculate_coefficients(angles,no_rho);
%         
%         % differentiation matrix
%         diff = sparse(1:nt,1:nt,-1*ones(1,nt),nt,nt+1) + sparse(1:nt,2:nt+1,1*ones(1,nt),nt,nt+1);
%         
%         G = build_operator_g(wavelets,diff,c1,c2,c3);
%         
%         
%         

% end
            
end