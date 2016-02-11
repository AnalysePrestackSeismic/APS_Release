function[]=run_batch_2d_algo_restart(batch_run_info_path,varargin)
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
%  
% if you dont want to scan, just want to run an algo in loop
% input_filepath = path of batch run info file
% varargin :  algorithm name, slurm_part, n_cores, then all algorithm parameters as strings

% Note: you dont need to give the wkey since these information will be
% extracted form the batch file tun info mat file. Also the wkey will be
% different for different lines.

%%
batch_run_info_path = batch_run_info_path{1};
batch_run_file_info=load(batch_run_info_path);       % Load batch run info path
rg=1;                                               % Initialize varargin counter to 0
algorithm_name=varargin{rg};rg=rg+1;                %get various parameters from varargin
slurm_part=varargin{rg};rg=rg+1;
n_cores=varargin{rg};rg=rg+1;
param_script=[varargin{rg}];rg=rg+1;                % Start making the tail script for parameters
for kk=rg:length(varargin)
    param_script=[param_script,' ',varargin{kk},];  % construct parmaeter script
end


for g=1:batch_run_file_info.nfiles
    id=strfind(batch_run_file_info.algo_run_qc{g},'wckey=');
    current_wckey=batch_run_file_info.algo_run_qc{g}((id+6):end);
    batch_run_file_info.wkey_orig{g}=current_wckey;
    text_script = ['/apps/gsc/matlab-library/development/maps/node/run_node_slurm_submit_neo_restart.sh /apps/matlab/v2013a/',' ',algorithm_name,' ',batch_run_file_info.job_meta_files{g},' ',slurm_part,' ',n_cores,' ',current_wckey,' ',param_script];
    [~,text_return]=system(text_script);
    text_return2=deblank(regexp(text_return,'\n','split'));
    batch_run_file_info.algo_run_qc_restart{g}=text_return2{1,9};
    pause(5);
    fprintf('submitted line : %g \n',g)
    
    
    
end

save(batch_run_info_path,'-struct','batch_run_file_info','-v7.3'); % Saves Seismic structure to mat file
disp(batch_run_info_path);


end