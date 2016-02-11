function[]=run_batch_2d_check_status(batch_run_info_path)
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
% input: batch_run_info_path = path to batch run info mat file
% prints status of each job
%%

altnodes = {'sblxpgcn192 ','sblxpghn100 ','sblxpgcn182 '};
fprintf('testing to contact head node or try alternate\n');

for hnloop = 1:1:size(altnodes,2)
    fprintf('testing to contact node %s \n',altnodes{hnloop});
    [sshstat,sshreply] = system(['ssh -x ',altnodes{hnloop},' date']);
    
    if sshstat == 0
        head_node = altnodes{hnloop};
        fprintf('slected node %s as head node\n',altnodes{hnloop});
        break;
    end
end

batch_run_file_info=load(batch_run_info_path);


% Display all the slurm commands fopr getting status reports
for p=1:batch_run_file_info.nfiles
    
    disp(batch_run_file_info.algo_run_qc{p});
    
end
disp('################################ STATUS REPORT #######################################')
%----------------- Processing the status of various jobs-----------
% batch_run_file_info.run_complete=ones(batch_run_file_info.nfiles,1);
for g=1:batch_run_file_info.nfiles
    job_meta=load(batch_run_file_info.job_meta_files{g});
    jobstatus = ones(size(job_meta.block_keys,1),1);
%     comm = [batch_run_file_info.algo_run_qc{g},' ','--state=FAILED'];
    comm = [batch_run_file_info.algo_run_qc{g}];
%     system(['ssh -x ',head_node,' ',comm]);
    [~,jobs] = system(['ssh -x ',head_node,' ',comm]);              % sssh into a node in the cluster and execute the sacct command
    jobcell =  strjust(regexp(jobs,'\n','split'),'left');           % process the text returned from the sacct command
    jobcell = deblank(jobcell); 
    if  size(jobcell,2) >7 
%         batch_run_file_info.run_complete(g) = 0;
        disp(batch_run_file_info.files_in.names{g});
        disp(jobcell{5});
        disp(jobcell{6});
        for k=7:size(jobcell,2)
            disp(jobcell{k});
        end
    end
    
      
end

% disp('----------List of Fully Complete Lines---------------')
% for g=1:batch_run_file_info.nfiles
%     
%     if batch_run_file_info.run_complete(g)==1
%         disp(batch_run_file_info.files_in.names{g});
%     end
% end


save(batch_run_info_path,'-struct','batch_run_file_info','-v7.3');
    
end