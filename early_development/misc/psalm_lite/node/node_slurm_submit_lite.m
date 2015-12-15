function [] = node_slurm_submit_lite(algorithm_name,i_block,n_blocks,seismic_mat_path,slurm_part,nodelist,number_cores,varargin)
% submit compiled matlab function slurm
% arguments cell array needs to be in correct order
% water_bottom_flatten(block_mat,process_files_mat,join_block_id)

    run_script_path = '/apps/gsc/matlab-mcode-beta/eslib/psalm_lite/algorithms/run_function.sh';
    matlab_path = '/apps/matlab/v2011a';
    
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

     for ii_block = i_block:1:n_blocks
        % make unique job directory in /localcache/mcr_cache_umask_friendly/ by submission through slurm 
%         if size(nodelist,2) > 0
%             for i_node = 1:1:size(nodelist,2)
%                 nodes
%             end                
%             slurm_job = sprintf('srun -p %s -c %s -J %s_block_%d %s %s %s %s %s %s %s %s',slurm_part,number_cores,algorithm_name,i_block,run_script_path,matlab_path,algorithm_name,seismic_mat_path,num2str(i_block),num2str(n_blocks),arguments,'&');
%             system(slurm_job);            
%         else    
            slurm_job = sprintf('srun -p %s -c %s -J %s_block_%d %s %s %s %s %s %s %s %s',slurm_part,number_cores,algorithm_name,ii_block,run_script_path,matlab_path,algorithm_name,seismic_mat_path,num2str(ii_block),num2str(n_blocks),arguments,'&');
            if ii_block == 1
                fprintf('%s\n',slurm_job);
            end
            system(slurm_job);
%         end
        fprintf('Function %s, job %d submitted to SLURM\n',algorithm_name,ii_block);
     end
end