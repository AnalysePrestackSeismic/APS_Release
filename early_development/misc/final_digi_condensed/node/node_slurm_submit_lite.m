function [] = node_slurm_submit_lite(algorithm_name,job_meta_path,slurm_part,n_cores,varargin)
%function [] = node_slurm_submit_lite(algorithm_name,n_blocks,seismic_mat_path,slurm_part,number_cores,varargin)
% submit compiled matlab function slurm
% arguments cell array needs to be in correct order
% water_bottom_flatten(block_mat,process_files_mat,join_block_id)

run_script_path = '/apps/gsc/matlab-mcode-beta/eslib/final_digi_condensed/algorithms/';
matlab_path = '/apps/matlab/v2011a';
job_meta = load(job_meta_path);
n_blocks = num2str(job_meta.n_blocks);
ncpus_per_node = 12;
% nodelist {tvlxpgcn70,tvlxpgcn71}

% function should be compiled as single threaded so a simple call to srun
% can be used
% Thames [01-16]
% UK1 [51-56, 58-98]
% All [01-16,51-56,58-98]
if strcmp(slurm_part,'Thames')
    n_nodes = 16;
    for i_node = 1:1:n_nodes
        slurm_purge = sprintf('srun -p %s -n 1 --exclusive %s %s','Thames','/apps/gsc/matlab-mcode-beta/eslib/psalm_lite/algorithms/purge_function.sh','&');
        system(slurm_purge);
    end
    pause(1);
elseif strcmp(slurm_part,'UK1')
    n_nodes = 48;
    for i_node = 1:1:n_nodes
        slurm_purge = sprintf('srun -p %s -n 1 --exclusive %s %s','UK1','/apps/gsc/matlab-mcode-beta/eslib/psalm_lite/algorithms/purge_function.sh','&');
        system(slurm_purge);
    end
    pause(1);
end

% Extra Arguments
arguments = '';
for i_arg = 1:1:length(varargin)
    if ~ischar(varargin{i_arg})
        varargin{i_arg} = num2str(varargin{i_arg});
    end
    arguments = sprintf('%s %s',arguments,char(varargin{i_arg}));
end

random_string = randi(1000000,1);
batch_script_path = [job_meta.output_dir,'script_job_',num2str(random_string)];
fid_batch = fopen(batch_script_path, 'w');
% fprintf(fid_batch, '%s\n', '#!/bin/sh');
% fprintf(fid_batch, '%s\n', ['#SBATCH --job-name=',algorithm_name]); 
% fprintf(fid_batch, '%s\n', ['#SBATCH -p ',slurm_part]);
% fprintf(fid_batch, '%s\n', '#SBATCH --error=job.%J.err');
% fprintf(fid_batch, '%s\n', '#SBATCH --output=job.%J.out');
% fprintf(fid_batch, '%s\n', ['#SBATCH --nodes=',num2str(n_nodes)]); 
% fprintf(fid_batch, '%s\n', ['#SBATCH --ntasks=',n_tasks]); 
% %fprintf(fid_batch, '%s\n', ['#SBATCH --cpus-per-task=1']);
% %fprintf(fid_batch, '%s\n', ['#SBATCH --ntasks-per-node=10']);
% fprintf(fid_batch, '%s\n', '#SBATCH --exclusive');
fprintf(fid_batch, 'umask 002\n');
fprintf(fid_batch, 'MCR_CACHE_ROOT=/localcache/mcr_cache_umask_friendly/${LOGNAME}_${RANDOM};\n');
fprintf(fid_batch, 'MCRROOT=/apps/matlab/v2011a/ ;\n');
fprintf(fid_batch, 'export MCR_CACHE_ROOT ;\n');
fprintf(fid_batch, 'export MCRROOT ;\n');
fprintf(fid_batch, 'mkdir -p ${MCR_CACHE_ROOT} ;\n');
fprintf(fid_batch, 'LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;\n');
fprintf(fid_batch, 'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;\n');
fprintf(fid_batch, 'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;\n');
fprintf(fid_batch, 'MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64 ;\n');
fprintf(fid_batch, 'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads ; \n');
fprintf(fid_batch, 'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server ;\n');
fprintf(fid_batch, 'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/client ;\n');
fprintf(fid_batch, 'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE} ;  \n');
fprintf(fid_batch, 'XAPPLRESDIR=${MCRROOT}/X11/app-defaults ;\n');
fprintf(fid_batch, 'export LD_LIBRARY_PATH;\n');
fprintf(fid_batch, 'export XAPPLRESDIR;\n');

for i_block = 1:1:str2double(n_blocks)
    % make unique job directory in /localcache/mcr_cache_umask_friendly/ by submission through slurm
    %         if size(nodelist,2) > 0
    %             for i_node = 1:1:size(nodelist,2)
    %                 nodes
    %             end
    %             slurm_job = sprintf('srun -p %s -c %s -J %s_block_%d %s %s %s %s %s %s %s %s',slurm_part,number_cores,algorithm_name,i_block,run_script_path,matlab_path,algorithm_name,seismic_mat_path,num2str(i_block),num2str(n_blocks),arguments,'&');
    %             system(slurm_job);
    %         else
    %slurm_job = sprintf('srun -p %s -c %s -J %s_block_%d %s %s %s %s %s %s %s %s',slurm_part,number_cores,algorithm_name,i_block,run_script_path,matlab_path,algorithm_name,seismic_mat_path,num2str(i_block),num2str(n_blocks),arguments,'&');
    % int_grad_inv_proj(job_meta_path,i_block)
    %job_script_path = [job_meta.output_dir,'script_job_',num2str(i_block),'.batch'];
    %fid_script = fopen(job_script_path, 'w');
    
    %slurm_job = [run_script_path,' ',algorithm_name,' ',job_meta_path,' ',num2str(i_block)];
    %fprintf(fid_script, '%s\n', slurm_job);
    fprintf(fid_batch, ['srun -p ',slurm_part,' -c ',n_cores,' -J ',algorithm_name,'_',num2str(i_block),' -o ',job_meta.output_dir,algorithm_name,'_',num2str(i_block),'.out ',run_script_path,algorithm_name,'/',algorithm_name,' ',job_meta_path,' ',num2str(i_block)]);
    %if i_block < str2double(n_blocks)
        fprintf(fid_batch, '%s\n', ' &');  
    %else
    %    fprintf(fid_batch, '\n');  
    %end
    %fclose(fid_script);
    %slurm_job = sprintf('srun -p %s -c %s -J %s_block_%d %s %s %s %s %s %s %s %s %s %s',slurm_part,number_cores,...
    %    algorithm_name,run_script_path,job_meta_path,i_block,matlab_path,algorithm_name,num2str(i_block),num2str(n_blocks),'100','&');
    %seismic_anomaly_spotter(job_meta_path,slice_file_rel_path,i_vol,i_slab,n_slab,window_length)
    %slurm_job = sprintf('srun -p %s -c %s -w %s -J %s_block_%d %s %s %s %s %s %s %s %s',slurm_part,number_cores,nodelist,algorithm_name,i_block,run_script_path,matlab_path,algorithm_name,seismic_mat_path,num2str(i_block),num2str(n_blocks),arguments,'&');
    
    % system(slurm_job);
    %         end
    % fprintf('Function %s, job %d submitted to SLURM\n',algorithm_name,i_block);
end

%fprintf(fid, '%s\n', ['echo "',n_blocks,' jobs submitted to SLURM partition ',slurm_part,'"']);
fprintf(fid_batch, '%s\n', 'exit');
fclose(fid_batch);

system(['chmod 777 ',batch_script_path]);
system(['ssh tvlxpgcn01 ',batch_script_path]);
end