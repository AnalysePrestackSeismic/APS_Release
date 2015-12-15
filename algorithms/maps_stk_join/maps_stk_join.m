function [] = maps_stk_join(job_meta_path,varargin)
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
%% Function Description
%  This function works on migrated gathers to figure out the residual
%  statics through a sophisticated trim statics approach . It can flatten
%  the gather and also produce a gather of the prestack shifts required. It
%  can also make a measurement of average absolute value of shift required
%  across the gather to flatten the event which we can the MICA (Migrated
%  Image confidence Attrbute). The process has higher resolution than
%  conventional trim statics approac and more robust to the scale inssues
%  associated with trim statics
%
% #####################Input:##############################
%   job_meta_path:      Path to job meta file
%
%   vargin:             Optional extra flags for zmin and zmax to limit the
%   trace length
% ######################Output:###########################
%   Void Output
%
% ###################### Writes to Disk:#######################
%   Everything gets written in directory : .../trimout/*. 
%   The  volumes written out items are:
%
% Can also output depending on outputflatgathers
%       Flattened Gathers (after applying shifts)
%
% sort out input parameters
%
% pick up argements to truncate the z range
ztrunc = 0;
lpi = 1;
if ~isempty(varargin)
    % test to make sure zmin and zmax were set
    if length(varargin) == 2
        zmin = str2double(varargin{1});
        zmax = str2double(varargin{2});
        ztrunc = 1;
    else
        error('must specify zmin and zmax on command line if using z limits');
    end
end

%
% if isempty(regexp(zsmooth,'il','once')) == 0
%     useselectemode = 1;
%     requiredinline =  str2double(regexprep(zsmooth,'il',''));
%     %requiredinline =  str2double(strrep(tottracerun,'il',''));
%     zsmooth = 20;
% else
%     zsmooth = str2double(zsmooth);
% end

% ============================================================================================================
% load the job meta file containing more parameters
job_meta = load(job_meta_path);
orig_sample_rate = (job_meta.s_rate/1000);


% add the history of jobs run and this one to the curent ebcdic
ebdichdr = ['merge '];
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
%
% for anything other than 2d this needs to have ratnos code to loop round the blocks in inline chunks
% and then an array index that makes the output linear in inlines not in
% crosslines and then an outer loop to read the next set of blocks in
% ============================================================================================================
% Read the data for this block
if isfield(job_meta, 'liveblocks')
    loopfin = size(job_meta.liveblocks,1);
    while lpi <= loopfin
        %for i_block = 1:1:str2double(n_blocks)
        
        i_block = num2str(job_meta.liveblocks(lpi));
        % Read the data for this block
        if lpi == 1
            if ztrunc == 1
                [~, results_out{2,2}, results_out{1,2}{1,1}, offset_read] = node_segy_read(job_meta_path,'1',i_block,zmin,zmax);
            else
                [~, results_out{2,2}, results_out{1,2}{1,1}, offset_read] = node_segy_read(job_meta_path,'1',i_block);
            end
        else
            if ztrunc == 1
                [~, tmp_results_out1, tmp_results_out2, tmp_offset_read] = node_segy_read(job_meta_path,'1',i_block,zmin,zmax);
            else
                [~, tmp_results_out1, tmp_results_out2, tmp_offset_read] = node_segy_read(job_meta_path,'1',i_block);
            end
            results_out{2,2} = [results_out{2,2} tmp_results_out1];
            results_out{1,2}{1,1} = [results_out{1,2}{1,1}; tmp_results_out2];
            offset_read = [offset_read tmp_offset_read];
        end
        
        % find unique offsets
        %offset = unique(offset_read);
        lpi = lpi + 1;
    end
else
    error('live blocks not defined in job_meta file');
end
%results_out{2,2} = reshape(results_out{2,2},size(results_out{2,2},1),length(offset),[]);


% check to make sure it read something if not exit
if isempty(results_out{2,2}) == 1 &&  isempty(results_out{1,2}{1,1}) == 1 %&& isempty(offset_read) == 1
    return
end

% drop data as required
% order is samples, angles, skeys
% results_out{2,2} = results_out{2,2}(1:maxzout,:,:);
% job_meta.n_samples{1} = maxzout;

[n_samples,n_traces] = size(results_out{2,2});
in_n_samples = n_samples;
%
% ============================================================================================================
% Make results meta information
%
output_dir = job_meta.output_dir;
% check to see if the directory exists
if exist(output_dir,'dir') == 0
    mkdir(output_dir);
else
    system('sleep 2');
    if exist(output_dir,'dir') == 0
        mkdir(output_dir);
    end
end
%
if job_meta.is_gather == 1
    % for the pre-stack dataset
    results_out{1,1} = 'Meta data for output files';
    %results_out{resultno,2}{1,1} = ilxl_read;
    results_out{1,2}{2,1} = offset_read';
    ebcstrtowrite = sprintf('%-3200.3200s',[results_out{1,1} '  ' ebdichdr '  ' tmpebc]);
    results_out{1,1} = ebcstrtowrite;
    results_out{1,3} = 'is_gather'; % 1 is yes, 0 is no
    
    % prestack output
    filename = regexprep( job_meta.files{1},'_block_[0-9]+.mat_orig_lite','');
    results_out{2,1} = strcat(filename,'_merge');
    %results_out{2,2} = zeros(n_samples,n_traces_gather,n_traces,'single');
    %commented out as written to above
    results_out{2,3} = 1;
    
else
    
    % for the stack output
    results_out{1,1} = 'Meta data for output files';
    %results_out{1,2}{1,1} = results_out{1,2}{1,1}(1:n_traces_gather:end,:);
    results_out{1,2}{2,1} = int32(zeros(n_traces,1));
    ebcstrtowrite = sprintf('%-3200.3200s',[results_out{1,1} '  ' ebdichdr '  ' tmpebc]);
    results_out{1,1} = ebcstrtowrite;
    results_out{1,3} = 'is_gather'; % 1 is yes, 0 is no
    %
    % stack result
    filename = regexprep( job_meta.files{1},'_block_[0-9]+.mat_orig_lite','');
    results_out{2,1} = strcat(filename,'_merge');
    results_out{2,2} = zeros(in_n_samples,n_traces,'single');
    results_out{2,3} = 0;
end
%
% ========================================================================
%
% % need to reshape the 2 3d matricies into 2d ones as the segy write just wants samples * total traces
% results_out{2,2} = reshape(results_out{2,2},size(results_out{2,2},1),length(offset),[]);
% results_out{2,2} = reshape(results_out{2,2},in_n_samples,[]);
% 
% if outputgathershifts == 1;
%     results_out{3,2} = reshape(results_out{3,2},in_n_samples,[]);
% end

i_block = str2double(i_block);
% write the stack dataset

if ztrunc == 1
    node_segy_write(results_out,i_block, orig_sample_rate, output_dir, zmin)
else
    node_segy_write(results_out,i_block, orig_sample_rate, output_dir)
end

file_name = fullfile(strcat(output_dir,results_out{2,1},'/'), sprintf('%s_block_%d%s',results_out{2,1},i_block,'.segy'));
fprintf('wrote file\n%s\n',file_name);

end