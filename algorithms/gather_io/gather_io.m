function [] = gather_io(job_meta_path,i_block,encrypt)
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
if encrypt == '1';
    %pw = '3230313542474f6374415753';
    % used for URY pw = '3230313542474f6374415753';
end
%%

job_meta = load(job_meta_path);
orig_sample_rate = (job_meta.s_rate/1000);

% add the history of jobs run and this one to the curent ebcdic
ebdichdr = ['trim statics '];
if isfield(job_meta,'comm_history')
    ebdichdr2 = job_meta.comm_history;
    tmpebc = ebdichdr2{size(ebdichdr2,1),2};
else
    ebdichdr2{1,2} = '';
    tmpebc = '';
end

for ebcii = (size(ebdichdr2,1)-1):-1:1
    tmpebcc = regexp(ebdichdr2{ebcii,2},'/','split');
    tmpebc = [tmpebc tmpebcc{1}  tmpebcc{end}];
end
tmpebc = sprintf('%-3200.3200s',tmpebc);
clear tmpebcc ebdichdr2;

% Read the data for this block
[~, results_out{2,2}, results_out{1,2}{1,1}, offset_read] = node_segy_read(job_meta_path,'1',i_block);

% should we update this for angle stacks?

offset = unique(offset_read);

% check to make sure it read something if not exit
if isempty(results_out{2,2}) == 1 &&  isempty(results_out{1,2}{1,1}) == 1 && isempty(offset_read) == 1
    return
end

[n_samples,n_traces_gather,n_traces] = size(results_out{2,2});
in_n_samples = n_samples;

% for the pre-stack dataset
results_out{1,1} = 'Meta data for output files';
%results_out{resultno,2}{1,1} = ilxl_read;
results_out{1,2}{2,1} = offset_read';
ebcstrtowrite = sprintf('%-3200.3200s',[results_out{1,1} '  ' ebdichdr '  ' tmpebc]);
results_out{1,1} = ebcstrtowrite;
results_out{1,3} = 'is_gather'; % 1 is yes, 0 is no

output_dir = [job_meta.output_dir,'gather_io/'];
% check to see if the directory exists
if exist(output_dir,'dir') == 0
    mkdir(output_dir);
end

% prestack output
filename = 'gaths';
results_out{2,1} = strcat(filename,'_encrypt_data');

%results_out{2,2} = zeros(n_samples,n_traces_gather,n_traces,'single');
%commented out as written to above
results_out{2,3} = 1;

i_block = str2double(i_block);
% write the pre-stack dataset

filename_out = node_segy_write(results_out,i_block, orig_sample_rate, output_dir);

%% Gather encryption
if encrypt == '1';    
    gpg_command = ['gpg --yes --batch --passphrase=',pw,' -c ',filename_out];
    [status,cmdout] = system(gpg_command);
    if status == 0
        rm_file = ['rm ',filename_out];
        system(rm_file);
        fprintf('Completed encryption with status \n');
    end    
end

%%

end