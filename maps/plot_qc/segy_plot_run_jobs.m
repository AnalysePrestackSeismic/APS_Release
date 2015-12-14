function segy_plot_run_jobs(job_meta_path,varargin)
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
%SEGY_PLOT_RUN_JOBS : Shows the status of the job-run on various blocks.
% Arguments: job_meta_path: Path of the job meta file
% optional argument wckey to select a different wckey from slurm
% e.g. segy_plot_run_jobs('/data/URY/segy/2013_pgs_uruguay_processing/full_area_final_deliverables_phase1_and_2/bg_matlab_ouput/job_meta/job_meta_29Oct2014.mat','jonesce_73590185320_trim_calculation')
%--------------------------------------------------------------------------
altnodes = {'sblxpgcn192 ','sblxpghn100 ','sblxpgcn182 '}; 
 
usessh = 0;
usecurwckey = 1;
job_meta = load(job_meta_path);                 % Load the job meta file


for i_arg = 1:(length(varargin))                                        % use all the arguments in varargin except the horizon path
    if ~ischar(varargin{i_arg})
        varargin{i_arg} = num2str(varargin{i_arg});                     %if varargin element not a character convert the numbers into string
    end
end

if length(varargin) < 1
    %get current slurm job
    current_wckey = job_meta.comm_history{end,1};   % Finds wkey associated with the last command in saved command history
else
    if strcmp(varargin{1},'neo')
        usessh = 1;
        current_wckey = job_meta.comm_history{end,1};   % Finds wkey associated with the last command in saved command history
        % try to see if the head node is contactable otherwise try
        % alternate nodes
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
        
        if  length(varargin) == 2
            current_wckey = varargin{2};
            usecurwckey = 0;
        end    
    elseif  length(varargin) == 1
        current_wckey = varargin{1};
        usecurwckey = 0;
    end
end
disp('current wckey is:');
disp(current_wckey);
datemonprev = datestr((now-100),29);            % Get the date 100 days ago 

scrsz = get(0,'ScreenSize');
figure('OuterPosition',[scrsz(3)/10 scrsz(4)/40 scrsz(3)/1.5 scrsz(4)/1.06]);
%---------------------------------------------------------------------------


% ----------------Draw all the blocks out---------------------------------

for i_block = 1:1:size(job_meta.block_keys,1)
cjxdata(:,i_block) = [job_meta.block_keys(i_block,1); job_meta.block_keys(i_block,2); job_meta.block_keys(i_block,2); job_meta.block_keys(i_block,1)];
cjydata(:,i_block) = [job_meta.block_keys(i_block,3); job_meta.block_keys(i_block,3); job_meta.block_keys(i_block,4); job_meta.block_keys(i_block,4)];

cdata(1,i_block,1) = 0.93;
cdata(1,i_block,2) = 0.93;
cdata(1,i_block,3) = 0.93;

end
%-------------------------------------------------------------------------


%zdata = ones(4,size(job_meta.block_keys,1));
%patch(cjxdata,cjydata,zdata,'w');
p = patch(cjxdata,cjydata,'w');
set(p,'EdgeColor',[0.93 0.93 0.93]);
set(p,'MarkerEdgeColor','none');
%set(p,'FaceColor','flat','CData',cdata)

%declare array for jobstatus
jobstatus = zeros(size(job_meta.block_keys,1),1);

jobopts = {'R',1;'CD',2;'F,NF,CA',3;'PD',4;'S',5};

%loop round getting all the jobstatus , did not get all at the same time to
%reduce length of string passed back from system
formatin = 'yyyy-mm-ddTHH:MM:SS';

if usecurwckey == 1
    
    for ii= size(jobopts,1):-1:1;
        
        jobcom = ['sacct --format=JobName%-40 --noheader --state=',jobopts{ii,1},' --allusers --starttime=',datemonprev,' --wckey=',current_wckey];
        fprintf('running command: %s\n',jobcom);
        
        if(usessh == 0)
            [~,jobs] = system(jobcom);
            fslrw = 1;
        else
            [~,jobs] = system(['ssh -x ',head_node,' ',jobcom]);
            fslrw = 5;
        end
        
        jobcell = deblank(regexp(jobs,'\n','split'));
        
        for cjii = fslrw:1:(size(jobcell,2)-1)
            tmpcell = regexp(jobcell{1,cjii}, '_','split');
            if ~isempty(regexpi(tmpcell{1,end},'[0-9]+'))
                jobstatus(str2double(tmpcell{1,end})) = jobopts{ii,2};
            end
        end
        
    end
    
else
    %declare array for jobstatus
    jobstatustmp = zeros(size(job_meta.block_keys,1),2);
    for ii= size(jobopts,1):-1:1;
        
        jobcom = ['sacct --format=Start,JobName%-40 --noheader --state=',jobopts{ii,1},' --allusers --starttime=',datemonprev,' --wckey=',current_wckey];
        fprintf('running command: %s\n',jobcom);
        if(usessh == 0)
            [~,jobs] = system(jobcom);
            fslrw = 1;
        else
            [~,jobs] = system(['ssh -x ',head_node,' ',jobcom]);
            fslrw = 5;
        end        
        
        %jobcell = deblank(regexp(jobs,'\n','split'));
        %jobcell =  strjust(jobcell,'left');
        jobcell =  strjust(regexp(jobs,'\n','split'),'left');
        jobcell = deblank(jobcell);
        
        for cjii = fslrw:1:(size(jobcell,2)-1)
            tmpcell2 =  regexp(jobcell{1,cjii}, ' ','split');
            tmpcell = regexp(tmpcell2{1,2}, '_','split');
            if ~isempty(regexpi(tmpcell{1,end},'[0-9]+'))
                
                if jobstatustmp(str2double(tmpcell{1,end}),2) == 0
                    if strcmp(tmpcell2{1,1},'Unknown')
                        jobstatustmp(str2double(tmpcell{1,end}),2) = now;
                        jobstatustmp(str2double(tmpcell{1,end}),1) = jobopts{ii,2};
                    else
                        jobstatustmp(str2double(tmpcell{1,end}),2) = datenum(tmpcell2{1,1},formatin);
                        jobstatustmp(str2double(tmpcell{1,end}),1) = jobopts{ii,2};
                    end
                else
                    if datenum(tmpcell2{1,1},formatin) > jobstatustmp(str2double(tmpcell{1,end}),2)
                        jobstatustmp(str2double(tmpcell{1,end}),2) = datenum(tmpcell2{1,1},formatin);
                        jobstatustmp(str2double(tmpcell{1,end}),1) = jobopts{ii,2};
                    end    
                end
            end
        end
        
    end
    jobstatus = jobstatustmp(1:end,1);
end


%hold all;
%set(p,'FaceColor',[0 0.1 0]);
loopfin = size(job_meta.liveblocks,1);
lpi = 1;
while lpi <= loopfin
    i_block = job_meta.liveblocks(lpi);
    cdata(1,i_block,1) = 0.8;
    cdata(1,i_block,2) = 0.835;
    cdata(1,i_block,3) = 0.93;
    
    lpi = lpi + 1;
end
%

% colour the completed jobs
% {'Running',1;'Failed,NodeFailed',2;'Completed',3;'Pending',4;'Suspended',5};

for i_block = 1:1:(size(jobstatus,1))
    %
    switch jobstatus(i_block)
        case 1 % running
            cdata(1,i_block,1) = 0.99;
            cdata(1,i_block,2) = 0.97;
            cdata(1,i_block,3) = 0.6;
        case 2 %  completed
            cdata(1,i_block,1) = 0;
            cdata(1,i_block,2) = 0.95;
            cdata(1,i_block,3) = 0;         
        case 3 % Failed          
            cdata(1,i_block,1) = 1;
            cdata(1,i_block,2) = 0;
            cdata(1,i_block,3) = 0;               
        case 4 % pending
            cdata(1,i_block,1) = 1;
            cdata(1,i_block,2) = 0.8;
            cdata(1,i_block,3) = 0.6;
        case 5 % suspended
            cdata(1,i_block,1) = 0;
            cdata(1,i_block,2) = 0;
            cdata(1,i_block,3) = 0;            
    end
end

set(p,'FaceColor','flat','CData',cdata) % apply the colours to the faces

% add text to show the block

%for i_block = 1:20:size(job_meta.block_keys,1)
loopfin = size(job_meta.liveblocks,1);
lpi = 1;
while lpi <= loopfin    
    i_block = job_meta.liveblocks(lpi);
    figlabels = text(cjxdata(1,i_block),cjydata(2,i_block),num2str(i_block),'FontSize',6);
    lpi = lpi + 10;
end

title([tmpcell{1,1:(end-1)}, '   completed = green, failed = red, running = yellow, pending = orange' ]);   % add a title

% to add a legend correctly need to make sperate patches for eac h type of
% job

% jobleg = {'R';'F,NF';'CD';'PD';'S'};
% leg = legend(jobleg,'Location','EastOutside');
% set(leg,'Color',[0.43 0.93 0.63]);
% %leg = legend('completed',running);

axis tight;


end