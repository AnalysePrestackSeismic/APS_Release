function [] = node_slurm_submit_overflow_restart(algorithm_name,job_meta_path,slurm_part,n_cores,current_wckey,varargin)
%function [] = node_slurm_submit_lite(algorithm_name,n_blocks,seismic_mat_path,slurm_part,number_cores,varargin)
%% Function Description
% Inputs:
    % algorithm_name:
    % job_meta_path:
    % slurm_part:
    % n_cores:
    % varargin:

% Outputs: Void

% Writes to Disk:
    % Log Files
    
% Notes: You may want to QC the running of the job using segy_plot_run_jobs.m

%%

% submit compiled matlab function slurm
% arguments cell array needs to be in correct order
% water_bottom_flatten(block_mat,process_files_mat,join_block_id)


slurm_part = 'Blades_Test';
n_cores = '2';
run_script_path = '/apps/gsc/matlab-library/development/maps/algorithms/';
matlab_path = '/apps/matlab/v2011a';
job_meta = load(job_meta_path);             % load job meta file

%===== section for the restart to select the failed jobs ==========

alljobs = job_meta.liveblocks;

limitnodes = '';
%limitnodes = ' --nodelist=tvlxpgcn[71-73,75] ';


%get current slurm job
%current_wckey = job_meta.comm_history{end,1};

datemonprev = datestr((now-100),29);

%declare array for jobstatus
jobstatus = zeros(size(job_meta.liveblocks,1),1);

%loop round getting all the jobstatus , did not get all at the same time to
%reduce length of string passed back from system

%jobcom = ['sacct --format=JobName%-40 --noheader --state=completed  --starttime=',datemonprev,' --wckey=',current_wckey];
%jobcom = ['sacct --format=JobName%-40 --noheader --state=F,NF --allusers --starttime=',datemonprev,' --wckey=',current_wckey];
jobcom = ['sacct --format=JobName%-40 --noheader --state=F --allusers --starttime=',datemonprev,' --wckey=',current_wckey];

% 
% [~,jobs] = system(jobcom);
% jobcell = deblank(regexp(jobs,'\n','split'));
% 
% for cjii = 1:1:(size(jobcell,2)-1)
%     tmpcell = regexp(jobcell{1,cjii}, '_','split');
%     jobstatus(cjii) = str2double(tmpcell{1,end});
% end


[~,jobs] = system(jobcom);
jobcell = deblank(regexp(jobs,'\n','split'));
cjiic = 1;
for cjii = 1:1:(size(jobcell,2)-1)
    tmpcell = regexp(jobcell{1,cjii}, '_','split');
    if ~isempty(regexpi(tmpcell{1,end},'[0-9]+'))
        %jobstatus(str2double(tmpcell{1,end})) = jobopts{ii,2};
        jobstatus(cjiic) = str2double(tmpcell{1,end});
        cjiic = cjiic + 1;
    end
end


job_meta.liveblocks = alljobs(ismember(alljobs,jobstatus) == 1);

% ==================================================================
randdir = num2str(floor(now*100000));
datemonprev = datestr((now-100),29);
[~,usrname] = system('whoami');
randdir = strcat(usrname,'_',randdir);      % Make directory name based on the current time and the name of the user

% Make directory to save results
logoutdir = strcat(job_meta.output_dir,'out_log_files/');
if exist(logoutdir,'dir') == 0
    mkdir(logoutdir)
    system(['chmod 777 ',logoutdir]);
end

% assume all blocks are run
%n_blocks = num2str(job_meta.n_blocks);
% load the list of blocks to run

% nodelist {tvlxpgcn70,tvlxpgcn71}

% function should be compiled as single threaded so a simple call to srun
% can be used
% Thames [01-16]
% UK1 [51-56, 58-98]
% All [01-16,51-56,58-98]
%
if strcmp(slurm_part,'Thames')
    n_nodes = 16;
    head_node = 'tvlxpgcn01 ';
    for i_node = 1:1:n_nodes
        slurm_purge = sprintf('srun -p %s -n 1 --exclusive %s %s %s %s','Thames',run_script_path,'purge_function.sh',randdir,' &');
        system(slurm_purge);
    end
    pause(1);
elseif strcmp(slurm_part,'UK1')
    n_nodes = 48;
    head_node = 'tvlxpgcn53 ';
    for i_node = 1:1:n_nodes
        slurm_purge = sprintf('srun -p %s -n 1 --exclusive %s %s %s %s','UK1',run_script_path,'purge_function.sh',randdir,' &');
        system(slurm_purge);
    end
    pause(1);
elseif strcmp(slurm_part,'Blades_Test')
    n_nodes = 20;
    head_node = 'tvlxpbws51 ';
    for i_node = 1:1:n_nodes
        slurm_purge = sprintf('srun -p %s -n 1 --exclusive %s %s %s %s','Blades_Test',run_script_path,'purge_function.sh',randdir,' &');
        system(slurm_purge);
    end
    pause(1);    
end

%--------------PREPARING ARGUMENTS ------------------
% All functions have first two arguments as function_name(job_meta_path,i_block)
% Extra Arguments

if strcmp(algorithm_name,'seismic_anomaly_spotter')
    n_blocks = varargin{1};
    n_samps_slice = floor(job_meta.n_samples{str2double(varargin{2})}/str2double(n_blocks));
    start_block = 1:n_samps_slice:job_meta.n_samples{str2double(varargin{2})};
    end_block = start_block+n_samps_slice-1;
end

if strcmp(algorithm_name,'wavelet_estimation')
    varargin{length(varargin)+1} = num2str(job_meta.liveblocks(1));         % Find the first live block from job_meta file and append varargin arraay
end

% if you are running int_grad_inv_proj and user has supplied a maxzout
% horizon mask, edit varargin)

if (strcmp(algorithm_name,'int_grad_inv_proj') && ~isempty(strfind(varargin{5},'/')))  
    horizon_path = varargin{5};                                             % Path of the horizon
    flag_horizon_mask = 1;                                                  % Flag that horizon mask will be used
    varargin = varargin(1:4);                                               % delete the last maxzout argument from varargin
    arguments_crop = '';
    for i_arg = 1:(length(varargin))                                        % use all the arguments in varargin except the horizon path
        if ~ischar(varargin{i_arg})
            varargin{i_arg} = num2str(varargin{i_arg});                     %if varargin element not a character convert the numbers into string
        end
        arguments_crop = sprintf('%s %s',arguments_crop,char(varargin{i_arg}));  % Append arguments_crop string by new entry from varargin 
    end
else
    flag_horizon_mask = 0;                                                  % Flag that horizon mask won't be used
    arguments = '';
    for i_arg = 1:1:length(varargin)
        if ~ischar(varargin{i_arg})
            varargin{i_arg} = num2str(varargin{i_arg});                     %if varargin element not a character convert the numbers into string
        end
        arguments = sprintf('%s %s',arguments,char(varargin{i_arg}));       % Append arguments string by new entry from varargin 
    end
end
%-----------------MAKING THE SCRIPT FILE FOR SUBMISSION-------------------

%--------Constructing Script Portion required for all algorithm----
random_string = randi(1000000,1);
batch_script_path = [job_meta.output_dir,'script_job_',randdir,'_',num2str(random_string)];
fid_batch = fopen(batch_script_path, 'w');
fprintf(fid_batch, '#!/bin/sh\n');
fprintf(fid_batch, 'umask 002\n');
fprintf(fid_batch, 'MCR_CACHE_ROOT=/localdata/localcache/%s\n',randdir);
%fprintf(fid_batch, 'MCR_CACHE_ROOT=/localcache/Paradigm/jonesce/%s\n',randdir);
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
%livetestexi = exist('job_meta','var');
%if livetestexi == 1;
%-----------------------------------------------------------------------

%%--------Constructing Script Portion specific to the algorithm used----
%--------------FOR SEISMIC ANOMALY SPOTTER_----------------------------
if strcmp(algorithm_name,'seismic_anomaly_spotter')
    for i_block = 1:1:str2double(n_blocks)        
        fprintf(fid_batch, ['sbatch -p ',slurm_part,' -c ',n_cores,' -J ',algorithm_name,'_',...
            num2str(i_block),' --wckey=',current_wckey,' -o ',logoutdir,algorithm_name,'_',...
            num2str(i_block),'.out <<EOF\n']);
        
        fprintf(fid_batch, '#!/bin/sh\n');
        fprintf(fid_batch, ['srun ',run_script_path,algorithm_name,'/',...
            algorithm_name,' ',job_meta_path,' ']);
        arguments_anom = [num2str(varargin{2}) ' ' num2str(varargin{3}) ' ' num2str(varargin{4}) ' ' num2str(varargin{5}) ' ' num2str(start_block(i_block)) ' ' num2str(end_block(i_block))];
        fprintf(fid_batch, arguments_anom);
        fprintf(fid_batch, ' \nEOF\n');
        fprintf(fid_batch, '%s\n','sleep 0.05');
    end  
%------------FOR ALL OTHER ALGORITHMS -----------------------------------
else
    if isfield(job_meta, 'liveblocks')    
        loopfin = size(job_meta.liveblocks,1);
        lpi = 1;
       
        if (flag_horizon_mask)
            block_maxzout = make_horizon_mask( job_meta_path,horizon_path,'0','0'); % find out maxzout for block from horizon mask
        end
        
        while lpi <= loopfin
        %for i_block = 1:1:str2double(n_blocks)
        
            i_block = job_meta.liveblocks(lpi);
                        
            %---------PATCH FOR HORIZON MASK----------------------
            
            if (flag_horizon_mask)
                i_block_maxzout = block_maxzout(lpi);                   % max zout for the block from horizon mask
                arguments = sprintf('%s %s',arguments_crop,char(num2str(i_block_maxzout)));   % Append arguments string by new entry from varargin 
            end
            %------------------------------------------------------
            lpi = lpi + 1;
            fprintf(fid_batch, ['sbatch -p ',slurm_part,' -c ',n_cores,' -J ',algorithm_name,'_',...
                num2str(i_block),' --wckey=',current_wckey,' -o ',logoutdir,algorithm_name,'_',num2str(i_block),'.out <<EOF\n']);
            fprintf(fid_batch, '#!/bin/sh\n');    
            fprintf(fid_batch, ['srun ',run_script_path,algorithm_name,'/',...
                algorithm_name,' ',job_meta_path,' ',num2str(i_block)]);
            
            fprintf(fid_batch, arguments);  
            fprintf(fid_batch, ' \nEOF\n');  
            fprintf(fid_batch, '%s\n','sleep 0.05');
            
        end
    else
        for i_block = 1:1:str2double(n_blocks)        
            fprintf(fid_batch, ['sbatch -p ',slurm_part,' -c ',n_cores,' -J ',algorithm_name,'_',...
                num2str(i_block),' --wckey=',current_wckey,' -o ',logoutdir,algorithm_name,'_',...
                num2str(i_block),'.out <<EOF\n']);
            fprintf(fid_batch, '#!/bin/sh\n');    
            fprintf(fid_batch, ['srun ',run_script_path,algorithm_name,'/',...
                algorithm_name,' ',job_meta_path,' ',num2str(i_block)]);
            fprintf(fid_batch, arguments);  
            fprintf(fid_batch, ' \nEOF\n');  
            fprintf(fid_batch, '%s\n','sleep 0.05'); 
        end
    end
end


% if strcmp(algorithm_name,'seismic_anomaly_spotter')
%     for i_block = 1:1:str2double(n_blocks)        
%         fprintf(fid_batch, ['srun -p ',slurm_part,' -c ',n_cores,' -J ',algorithm_name,'_',...
%             num2str(i_block),' --wckey=',randdir,'_',algorithm_name,' -o ',job_meta.output_dir,algorithm_name,'_',...
%             num2str(i_block),'.out ',run_script_path,algorithm_name,'/',...
%             algorithm_name,' ',job_meta_path,' ']);
%         arguments_anom = [num2str(varargin{2}) ' ' num2str(varargin{3}) ' ' num2str(varargin{4}) ' ' num2str(start_block(i_block)) ' ' num2str(end_block(i_block))];
%         fprintf(fid_batch, arguments_anom);
%         fprintf(fid_batch, ' %s\n', ' &');
%         fprintf(fid_batch, '%s\n','sleep 0.05');
%     end   
% else
%     if isfield(job_meta, 'liveblocks')    
%         loopfin = size(job_meta.liveblocks,1);
%         lpi = 1;
%         %while lpi <= 20
%         while lpi <= loopfin
%         %for i_block = 1:1:str2double(n_blocks)
%             i_block = job_meta.liveblocks(lpi);
%             lpi = lpi + 1;
%             fprintf(fid_batch, ['srun -p ',slurm_part,' -c ',n_cores,' -J ',algorithm_name,'_',...
%                 num2str(i_block),' --wckey=',randdir,'_',algorithm_name,' -o ',job_meta.output_dir,algorithm_name,'_',...
%                 num2str(i_block),'.out ',run_script_path,algorithm_name,'/',...
%                 algorithm_name,' ',job_meta_path,' ',num2str(i_block)]);
%             fprintf(fid_batch, arguments);  
%             fprintf(fid_batch, ' %s\n', ' &');  
%             fprintf(fid_batch, '%s\n','sleep 0.05');        
% 
%         end
%     else
%         for i_block = 1:1:str2double(n_blocks)        
%             fprintf(fid_batch, ['srun -p ',slurm_part,' -c ',n_cores,' -J ',algorithm_name,'_',...
%                 num2str(i_block),' --wckey=',randdir,'_',algorithm_name,' -o ',job_meta.output_dir,algorithm_name,'_',...
%                 num2str(i_block),'.out ',run_script_path,algorithm_name,'/',...
%                 algorithm_name,' ',job_meta_path,' ',num2str(i_block)]);
%             fprintf(fid_batch, arguments);
%             fprintf(fid_batch, ' %s\n', ' &');
%             fprintf(fid_batch, '%s\n','sleep 0.05');
%         end
%     end
% end



%fprintf(fid, '%s\n', ['echo "',n_blocks,' jobs submitted to SLURM partition ',slurm_part,'"']);
fprintf(fid_batch, '%s\n', 'exit');
fclose(fid_batch);


% % Add processing information to job meta 
% if (flag_horizon_mask)
%    
%    arguments = sprintf('%s %s',arguments_crop,char(horizon_path));   % Append arguments string by horizon_path if horizon mask used.
% end
% 
% newcommand = {strcat(randdir,'_',algorithm_name),['node_slurm_submit ',algorithm_name,' ',job_meta_path,' ',slurm_part,' ',n_cores,' ',arguments]};
% %newcommand{end,1}
% if isfield(job_meta,'comm_history')
%     job_meta.comm_history = [job_meta.comm_history;newcommand];
% else
%     job_meta.comm_history = newcommand;
% end
% save(job_meta_path,'-struct','job_meta','-v7.3');


% Make the script file executable and log into node 1 to submit the job
system(['chmod 777 ',batch_script_path]);
system(['ssh ',head_node,batch_script_path,' &']);
system('exit');
fprintf('run the command below to see the job status\n');
fprintf('sacct --format=JobName%%-30,jobid,elapsed,Start,End,ncpus,ntasks,state,wckey%%-40,nodelist --allusers  --starttime=%s  --wckey=%s%s%s\n',datemonprev,randdir,'_',algorithm_name);

end