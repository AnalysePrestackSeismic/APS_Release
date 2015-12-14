function [] = node_slurm_submit_restart(algorithm_name,job_meta_path,slurm_part,n_cores,current_wckey,varargin)
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
% node_slurm_submit: function to submit jobs to cluster. A specific
% algorithm is specified and this each block is processed in parallel on a
% cluster
%   Arguments:
%       algorithm_name = the algorithm that should be submitted to the
%       cluster. This should have already been compiled using compile_function
%       job_meta_path = path to job_meta .mat file contain meta data
%       slurm_part = slurm partition to use
%       n_cores = specify how many cores the job requires. Used by -c flag
%       in slurm. Helpful for jobs that require more memory or that
%       multi-thread
%       current_wckey = the current job key from slurm, found in job.meta
%       or printed out by node_slurm_submit or by segy_plot_run_jobs
%       varargin = algorithm specific flags to control arguments. By
%       default all algorithms require the following arguments:
%           job_meta_path
%           i_block
%           custom arguments can be added '30','1','30','0','1'
%
%   Outputs:
%       none
%
%   Writes to Disk:
%       Log Files
% Notes: You may want to QC the running of the job using segy_plot_run_jobs.m
%%
%============================================================================================================================
% submit compiled matlab function to slurm
% arguments cell array needs to be in correct order
% constants
run_script_path = '/apps/gsc/matlab-library/development/maps/algorithms/';
matlab_path = '/apps/matlab/v2011a';
%
findfailed = 0;  % if this is set to anything other than 0 it finds only the failed jobs, otherwise finds the completed jobs and reruns anything that has not completed
%=============================================================================================================================
%%
job_meta = load(job_meta_path);

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
%
check_path = regexp(current_wckey,'/');
if ~isempty(check_path) && check_path(1) == 1 % load file
    job_meta.liveblocks = dlmread(current_wckey);
else
    if findfailed == 0
        jobcom = ['sacct --format=JobName%-40 --noheader --state=completed --allusers --starttime=',datemonprev,' --wckey=',current_wckey];
        %jobcom = ['sacct --format=JobName%-40 --noheader --state=completed --allusers --starttime=2014-12-29 --wckey=',current_wckey];
        
    else
        jobcom = ['sacct --format=JobName%-40 --noheader --state=F --allusers --starttime=',datemonprev,' --wckey=',current_wckey];
        %jobcom = ['sacct --format=JobName%-40 --noheader --state=F --allusers --starttime=2015-01-09 --wckey=',current_wckey];
    end
    [~,jobs] = system(jobcom);
    jobcell = deblank(regexp(jobs,'\n','split'));
    
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
    
    if findfailed == 0
        job_meta.liveblocks = alljobs(ismember(alljobs,jobstatus) == 0);
    else
        job_meta.liveblocks = alljobs(ismember(alljobs,jobstatus) == 1);
    end
end
% 
% [~,jobs] = system(jobcom);
% jobcell = deblank(regexp(jobs,'\n','split'));
% 
% for cjii = 1:1:(size(jobcell,2)-1)
%     tmpcell = regexp(jobcell{1,cjii}, '_','split');
%     jobstatus(cjii) = str2double(tmpcell{1,end});
% end



randdir = num2str(floor(now*100000));
datemonprev = datestr((now-100),29);
[~,usrname] = system('whoami');
randdir = strcat(usrname,'_',randdir);   % Make directory name based on the current time and the name of the user

logoutdir = strcat(job_meta.output_dir,'out_log_files/');

% assume all blocks are run
%n_blocks = num2str(job_meta.n_blocks);
% load the list of blocks to run

% nodelist {tvlxpgcn70,tvlxpgcn71}

% function should be compiled as single threaded so a simple call to srun
% can be used
% Thames [01-16]
% UK1 [51-56, 58-98]
% AllDC [01-16,51-56,58-98]
if strcmp(slurm_part,'Thames') || strcmp(slurm_part,'AllDC')
    n_nodes = 16;
    %head_node = 'tvlxpgcn01 ';
    head_node = 'tvlxpgcn100 ';
    for i_node = 1:1:n_nodes
        slurm_purge = sprintf('srun -p %s -n 1 --exclusive %s %s %s %s','Thames',run_script_path,'purge_function.sh',randdir,' &');
        system(slurm_purge);
    end
    pause(1);
elseif strcmp(slurm_part,'UK1')
    n_nodes = 48;
    %head_node = 'tvlxpgcn53 ';
    head_node = 'tvlxpgcn101 ';
    for i_node = 1:1:n_nodes
        slurm_purge = sprintf('srun -p %s -n 1 --exclusive %s %s %s %s','UK1',run_script_path,'purge_function.sh',randdir,' &');
        system(slurm_purge);
    end
    pause(1);
end
% All functions have first two arguments as function_name(job_meta_path,i_block)
% Extra Arguments

if strcmp(algorithm_name,'seismic_anomaly_spotter')
    n_blocks = varargin{1};
    n_samps_slice = floor(job_meta.n_samples{str2double(varargin{2})}/str2double(n_blocks));
    start_block = 1:n_samps_slice:job_meta.n_samples{str2double(varargin{2})};
    end_block = start_block+n_samps_slice-1;
end

if strcmp(algorithm_name,'wavelet_estimation')
    varargin{length(varargin)+1} = num2str(job_meta.liveblocks(1));        % Find the first live block from job_meta file and append varargin arraay
end

if strcmp(algorithm_name,'structural_tensor_dip')
    %scale_sigma = varargin{1};
    scale_sigma = varargin{3};
    sigma = varargin{2};
    aperture = num2str(str2num(scale_sigma)*str2num(sigma));
    %add_aperture_to_job_meta(job_meta_path,aperture);
    %node_slurm_submit('structural_tensor_dip','/data/TZA/segy/2013_kusini_inboard/pgs_enhanced_volume/structural_tensors_dip/job_meta/job_meta_07Oct2014.mat','UK1','4','1','3','1','3000')
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
fprintf(fid_batch, 'MCR_CACHE_ROOT=/localcache/mcr_cache_umask_friendly/%s\n',randdir);
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
        fprintf(fid_batch, ['sbatch -p ',slurm_part,limitnodes,' -c ',n_cores,' -J ',algorithm_name,'_',...
            num2str(i_block),' --wckey=',current_wckey,' -o ',logoutdir,algorithm_name,'_',...
            num2str(i_block),'.out <<EOF\n']);
        
        fprintf(fid_batch, '#!/bin/sh\n');
        fprintf(fid_batch, ['srun ',run_script_path,algorithm_name,'/',...
            algorithm_name,' ',job_meta_path,' ']);
        arguments_anom = [num2str(varargin{2}) ' ' num2str(varargin{3}) ' ' num2str(varargin{4}) ' ' num2str(start_block(i_block)) ' ' num2str(end_block(i_block))];
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
            fprintf(fid_batch, ['sbatch -p ',slurm_part,limitnodes,' -c ',n_cores,' -J ',algorithm_name,'_',...
                num2str(i_block),' --wckey=',current_wckey,' -o ',logoutdir,algorithm_name,'_',...
                num2str(i_block),'.out <<EOF\n']);
            fprintf(fid_batch, '#!/bin/sh\n');    
            fprintf(fid_batch, ['srun ',run_script_path,algorithm_name,'/',...
                algorithm_name,' ',job_meta_path,' ',num2str(i_block)]);
            fprintf(fid_batch, arguments);  
            fprintf(fid_batch, ' \nEOF\n');  
            fprintf(fid_batch, '%s\n','sleep 0.05');        
        end
    else
        for i_block = 1:1:str2double(n_blocks)        
            fprintf(fid_batch, ['sbatch -p ',slurm_part,limitnodes,' -c ',n_cores,' -J ',algorithm_name,'_',...
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
%             num2str(i_block),' --wckey=',current_wckey,' -o ',job_meta.output_dir,algorithm_name,'_',...
%             num2str(i_block),'.out ',run_script_path,algorithm_name,'/',...
%             algorithm_name,' ',job_meta_path,' ']);
%         arguments_anom = [num2str(varargin{2}) ' ' num2str(varargin{3}) ' ' num2str(varargin{4}) ' ' num2str(start_block(i_block)) ' ' num2str(end_block(i_block))];
%         fprintf(fid_batch, arguments_anom);
%         fprintf(fid_batch, ' %s\n', ' &');
%         fprintf(fid_batch, '%s\n','sleep 0.02');
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
%                 num2str(i_block),' --wckey=',current_wckey,' -o ',job_meta.output_dir,algorithm_name,'_',...
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
%                 num2str(i_block),' --wckey=',current_wckey,' -o ',job_meta.output_dir,algorithm_name,'_',...
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
fprintf('sacct --format=JobName%%-30,jobid,elapsed,Start,End,ncpus,ntasks,state,wckey%%-40,nodelist --allusers --starttime=%s  --wckey=%s%s%s\n',datemonprev,randdir,'_',algorithm_name);

end