function [] = node_slurm_submit_restart(algorithm_name,job_meta_path,slurm_part,n_cores,varargin)
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
% submit compiled matlab function slurm
% arguments cell array needs to be in correct order
% water_bottom_flatten(block_mat,process_files_mat,join_block_id)

run_script_path = '/apps/gsc/matlab-library/final_digi_condensed/algorithms/';
matlab_path = '/apps/matlab/v2011a';
job_meta = load(job_meta_path);

restartjobs = dlmread('/data/TZA/segy/2013_kusini_inboard/pgs_standard_volume/angle_stacks/final_digi_out/job_meta/restart_wavelet.jobs');
noofrestarts = size(restartjobs,1);

% assume all blocks are run
%n_blocks = num2str(job_meta.n_blocks);
% load the list of blocks to run

% nodelist {tvlxpgcn70,tvlxpgcn71}

% function should be compiled as single threaded so a simple call to srun
% can be used
% Thames [01-16]
% UK1 [51-56, 58-98]
% All [01-16,51-56,58-98]
if strcmp(slurm_part,'Thames')
    n_nodes = 16;
    head_node = 'tvlxpgcn01 ';
    for i_node = 1:1:n_nodes
        slurm_purge = sprintf('srun -p %s -n 1 --exclusive %s %s','Thames','/apps/gsc/matlab-mcode-beta/eslib/psalm_lite/algorithms/purge_function.sh','&');
        system(slurm_purge);
    end
    pause(1);
elseif strcmp(slurm_part,'UK1')
    n_nodes = 48;
    head_node = 'tvlxpgcn53 ';
    for i_node = 1:1:n_nodes
        slurm_purge = sprintf('srun -p %s -n 1 --exclusive %s %s','UK1','/apps/gsc/matlab-mcode-beta/eslib/psalm_lite/algorithms/purge_function.sh','&');
        system(slurm_purge);
    end
    pause(1);
end
% All functions have first two arguments as function_name(job_meta_path,i_block)
% Extra Arguments

if strcmp(algorithm_name,'seismic_anomaly_spotter')
    n_blocks = varargin{2};
end

if strcmp(algorithm_name,'wavelet_estimation')
    varargin{length(varargin)+1} = num2str(job_meta.liveblocks(1));
end

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

%if exist('job_meta.liveblocks','var');
    loopfin = noofrestarts;
    lpi = 1;
    %while lpi <= 20
    while lpi <= loopfin
    %for i_block = 1:1:str2double(n_blocks)
        i_block = job_meta.liveblocks(lpi);
        lpi = lpi + 1;
        fprintf(fid_batch, ['srun -p ',slurm_part,' -c ',n_cores,' -J ',algorithm_name,'_',...
            num2str(i_block),' -o ',job_meta.output_dir,algorithm_name,'_',...
            num2str(i_block),'.out ',run_script_path,algorithm_name,'/',...
            algorithm_name,' ',job_meta_path,' ',num2str(i_block)]);
        fprintf(fid_batch, arguments);  
        fprintf(fid_batch, ' %s\n', ' &');
        fprintf(fid_batch, '%s\n','sleep 0.01');
    end
% else
%     for i_block = 1:1:str2double(n_blocks)        
%         fprintf(fid_batch, ['srun -p ',slurm_part,' -c ',n_cores,' -J ',algorithm_name,'_',...
%             num2str(i_block),' -o ',job_meta.output_dir,algorithm_name,'_',...
%             num2str(i_block),'.out ',run_script_path,algorithm_name,'/',...
%             algorithm_name,' ',job_meta_path,' ',num2str(i_block)]);
%         fprintf(fid_batch, arguments);
%         fprintf(fid_batch, ' %s\n', ' &');
%     end
% end

%fprintf(fid, '%s\n', ['echo "',n_blocks,' jobs submitted to SLURM partition ',slurm_part,'"']);
fprintf(fid_batch, '%s\n', 'exit');
fclose(fid_batch);

% Make the script file executable and log into node 1 to submit the job
system(['chmod 777 ',batch_script_path]);
system(['ssh ',head_node,batch_script_path,' &']);
system('exit');
end