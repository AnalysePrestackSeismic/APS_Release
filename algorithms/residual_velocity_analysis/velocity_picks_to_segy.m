function [ output_args ] = velocity_picks_to_segy( job_meta, directory_path,file_stub,srate,outdir,output_type )
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
%VELOCITY_PICKS_TO_SEGY Read .mat files of RMS picks and make segy volume
%   job_meta:       meta data for gathers input to residual velocity analysis
%   directory path: location of .mat files with velocity picks
%   file_stub:      string to identify the .mat files to process
%   srate:          output sample rate
%   outdir:         directory for segy file. Output filenames will use
%                   file_stub plus _vrms/vint.segy
%   output_type:    vrms/vint/both

% convert srate number if supplied as string

if ischar(srate)
    srate = str2num(srate);
end

% add an underscore to the file_stub for output file if necessary

if strcmp(file_stub(end),'_')
    file_stub_out = file_stub;
else
    file_stub_out = [file_stub,'_'];
end

% add a forward slash to directory path if necessary

if ~strcmp(outdir(end),'/')
    outdir = [outdir,'/'];
end

vrmsfile = [outdir,file_stub_out,'vrms.segy'];
vintfile = [outdir,file_stub_out,'vint.segy'];

switch output_type
    case 'vrms'
        vrms_out = true;
        vint_out = false;
    case 'vint'
        vrms_out = false;
        vint_out = true;
    case 'both'
        vrms_out = true;
        vint_out = true;
end


% find inline/xline ranges from job_meta

gathers = load(job_meta);

tracelen = gathers.n_samples{1} * gathers.s_rate./1000;

n_vel_samples = ceil(tracelen/srate);

vel_tracelen = n_vel_samples * srate;



% find input files & scan to find inline/xline ranges

[scan_files nscanfiles] = directory_scan(cellstr(directory_path),file_stub);

% filter the file list for .mat files only & get list of blocks
% (there must be a neater way of getting the block number)

nfiles = 0;
live_blocks = zeros(nscanfiles,1);

for ff = 1:nscanfiles
    if strfind(scan_files.names{ff},'.mat')
        nfiles = nfiles + 1;
        input_files{nfiles} = scan_files.names{ff};
        tmp = regexp(input_files{nfiles},'block','split');
        tmp2 = regexp(tmp{2},'.mat','split');
        live_blocks(nfiles) = str2num(tmp2{1});
        clear tmp tmp2;
    end
end

live_blocks = live_blocks(1:nfiles,1);

live_block_keys = gathers.block_keys(live_blocks,:);
inl_ranges = unique(live_block_keys(:,1:2),'rows'); % find unique inline ranges



for rr = 1:size(inl_ranges,1)
    % find the files where the min and max inline is the same as the range
    % for this time round the loop
    file_idx = find((live_block_keys(:,1)-inl_ranges(rr,1)) + (live_block_keys(:,2)-inl_ranges(rr,2)) ==0);
    
    % min_xl/max_xl/num_xls/min_il/max_il/num_ils refer to inline/xline
    % ranges in the swath which share the same il range
    
    min_xl = min(live_block_keys(file_idx,3));
    max_xl = max(live_block_keys(file_idx,4));
    num_xls = 1+(max_xl-min_xl)/gathers.skey_inc;
    min_il = inl_ranges(rr,1);
    max_il = inl_ranges(rr,2);
    num_ils = 1+(max_il-min_il)/gathers.pkey_inc;
    vel_traces = zeros(n_vel_samples,num_xls*num_ils);
    
    for ff = 1:size(file_idx,1);
        vels = load([directory_path,'/',input_files{file_idx(ff)}]);
        
        % fix any zeros in the picks to be at the first sample
        % shouldn't be necessary now have changed the picking code
        
        vels.output_matrix(vels.output_matrix(:,3)==0,3) = gathers.s_rate./1000; 
               
        % interpolate onto regular z grid then reshape into a stream of
        % traces, keeping an index of inline/xline for each trace
        
        ilxl = unique(vels.output_matrix(:,1:2),'rows');
        ilxl_idx = 1+((ilxl(:,1)-min_il)./gathers.pkey_inc).*num_xls+(ilxl(:,2)-min_xl)./gathers.skey_inc;
        
       
        ntraces = size(ilxl,1);
        
        interp_times = repmat([srate:srate:vel_tracelen]',ntraces,1);
        
        interp_times_tmp = repmat(ilxl_idx,1,n_vel_samples)';
        
        interp_times = interp_times + interp_times_tmp(:).*vel_tracelen;
        
        vels.output_matrix(:,3) = vels.output_matrix(:,3) + (1 + (((vels.output_matrix(:,1)-min_il)./gathers.pkey_inc).*num_xls ... 
            + (vels.output_matrix(:,2)-min_xl)./gathers.skey_inc)) .* vel_tracelen;
        
        % check for increasing times
        [min_diff min_idx] = min(diff(vels.output_matrix(:,3)));
        
        if min_diff <=0 
            error(['Pick times not increasing in input file ',input_files{file_idx(ff)},' at line ',num2str(min_idx)]);
        end
        
        vel_traces_block = interp1(vels.output_matrix(:,3),vels.output_matrix(:,4),interp_times);
        vel_traces_block = reshape(vel_traces_block,n_vel_samples,ntraces);
        
        % put velocity traces into order
        
        vel_traces(:,ilxl_idx) = vel_traces_block;
         
        clear interp_times interp_times_tmp vel_traces_block ilxl ilxl_idx vels;
        
    end
        % write swath to segy
        
        % make array of inline/xline/offset for segy
        
        [ilgrid,xlgrid] = meshgrid(min_il:gathers.pkey_inc:max_il,min_xl:gathers.skey_inc:max_xl);
        headers = [ilgrid(:),xlgrid(:),zeros(num_ils*num_xls,1)]; 
    
    if rr==1
        mode = 'overwrite';
    else
        mode = 'append';
    end
    
    disp(['Writing inlines ',num2str(min_il),' to ',num2str(max_il),'.']);
    
    if vrms_out
        mat_to_segy(mode,'min',vel_traces,headers,srate*1000,'',vrmsfile);
    end
    
    if vint_out
        % make interval velocity volume
        vel_traces_vint = bsxfun(@vrms_to_vint,vel_traces,[srate:srate:vel_tracelen]');
        mat_to_segy(mode,'min',vel_traces_vint,headers,srate*1000,'',vintfile);
    end
end


    
    
end

    


