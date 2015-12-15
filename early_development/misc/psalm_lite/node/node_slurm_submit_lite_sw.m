function [] = node_slurm_submit_lite_sw(algorithm_name,n_blocks,seismic_mat_path,varargin)
% submit compiled matlab function slurm
% arguments cell array needs to be in correct order
% water_bottom_flatten(block_mat,process_files_mat,join_block_id)

run_script_path = '/apps/gsc/matlab-mcode-beta/eslib/psalm_lite/algorithms/run_function.sh';
matlab_path = '/apps/matlab/v2011a';

% function should be compiled as single threaded so a simple call to srun
% can be used
% Thames [01-16]
% UK1 [51-56, 58-98]
% All [01-16,51-56,58-98]
slurm_part = 'All';
n_nodes = 64;
for i_node = 1:1:n_nodes
    slurm_purge = sprintf('srun -p %s -n 1 --exclusive %s %s',slurm_part,'/apps/gsc/matlab-library/psalm_v0.1/algorithms/purge_function.sh','&');
    system(slurm_purge);
end
pause(5);

% Extra Arguments
arguments = '';
for i_arg = 1:1:length(varargin)
    if ~ischar(varargin{i_arg})
        varargin{i_arg} = num2str(varargin{i_arg});
    end
    arguments = sprintf('%s %s',arguments,char(varargin{i_arg}));
end


% failed_jobs = [1:12,97,101,141,142,182,183];
%for i_block = 1:1:n_blocks
for i_block = 1:1:n_blocks
    %for i_block = 62:1:73
    
    %    for ii = 1:1:size(failed_jobs(:))
    %        i_block = failed_jobs(ii)
    slurm_job = sprintf('srun -c 12 -p %s -J %s_block_%d %s %s %s %s %s %s %s %s',slurm_part,algorithm_name,i_block,run_script_path,matlab_path,algorithm_name,seismic_mat_path,num2str(i_block),num2str(n_blocks),arguments,'&');
    %slurm_job = sprintf('srun -c 12 -p %s -J %s_block_%d %s %s %s %s %s %s %s %s',slurm_part,algorithm_name,i_block,run_script_path,matlab_path,algorithm_name,seismic_mat_path,num2str(i_block),num2str(n_blocks),arguments,'&');
    %    slurm_job = sprintf('srun -x /apps/gsc/matlab-mcode-beta/eslib/psalm_lite/node/nodelist -c 1 -p %s -J %s_block_%d %s %s %s %s %s %s %s %s',slurm_part,algorithm_name,i_block,run_script_path,matlab_path,algorithm_name,seismic_mat_path,num2str(i_block),num2str(n_blocks),arguments,'&');
    
    system(slurm_job);
    fprintf('Function %s, job %d submitted to SLURM\n',algorithm_name,i_block);
end
end