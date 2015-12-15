function [] = node_slurm_submit(algorithm_name,n_blocks,path_for_blocks,varargin)
% submit compiled matlab function slurm
% arguments cell array needs to be in correct order
% water_bottom_flatten(block_mat,process_files_mat,join_block_id)

run_script_path = strcat('/apps/gsc/matlab-library/psalm_v0.1/algorithms/run_function.sh');
%func_path = strcat('/apps/gsc/matlab-library/psalm_v0.1/algorithms/',algorithm_name);
matlab_path = '/apps/matlab/v2011a';

% function should be compiled as single threaded so a simple call to srun
% can be used
% Thames [01-16]
% UK1 [51-56, 58-98]
% All [01-16,51-56,58-98]
% slurm_part = 'Thames';
slurm_part = 'UK1';
n_nodes = 48;
for i_node = 1:1:n_nodes
    slurm_purge = sprintf('srun -p %s -n 1 --exclusive %s %s',slurm_part,'/apps/gsc/matlab-library/psalm_v0.1/algorithms/purge_function.sh','&');
    system(slurm_purge);    
end
pause(5);
process_mat = strcat(path_for_blocks,algorithm_name,'_process_files.mat');

% Extra Arguments
arguments = '';
for i_arg = 1:1:length(varargin)
    if ~ischar(varargin{i_arg})
       varargin{i_arg} = num2str(varargin{i_arg});
    end        
    arguments = sprintf('%s %s',arguments,char(varargin{i_arg})); % --output /
end

for i_block = 1:1:n_blocks
    
    block_mat_all = strcat(path_for_blocks,algorithm_name,'_positions_block_all.mat');
    
    block_mat = strcat(path_for_blocks,algorithm_name,'_positions_block_',num2str(i_block),'.mat');
    
    %slurm_job = sprintf('srun -p %s -n 1 --exclusive -t %d %s %s %s %s %s %s %s %s',slurm_part,'0',run_script_path,matlab_path,algorithm_name,block_mat_all,block_mat,process_mat,arguments,'&');
    % slurm_job = sprintf('srun -c 12 -p %s %s %s %s %s %s %s %s %s',slurm_part,run_script_path,matlab_path,algorithm_name,block_mat_all,block_mat,process_mat,arguments,'&');
    slurm_job = sprintf('srun -c 12 -p %s %s %s %s %s %s %s %s %s',slurm_part,run_script_path,matlab_path,algorithm_name,block_mat_all,block_mat,process_mat,arguments,'&');
    % - c number of cpus that the job will use
    
    system(slurm_job);
    %pause(0.5)
    fprintf('Function %s, job %d submitted to SLURM\n',algorithm_name,i_block);

end

% Load sview

end