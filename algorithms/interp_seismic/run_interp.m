function run_interp(job_meta_path)
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

job_meta = load(job_meta_path);

    for i_block = 1:1:size(job_meta.liveblocks,1)
        interp_seismic(job_meta_path,num2str(job_meta.liveblocks(i_block)),'20','20');
    end

end

function [] = interp_seismic(job_meta_path,i_block,pkey_inc_factor,skey_inc_factor)

job_meta = load(job_meta_path);
pkey_inc_factor = str2double(pkey_inc_factor);
skey_inc_factor = str2double(skey_inc_factor);

% check inline crossline numbering
pkey_inc_mode = mode(job_meta.pkey_inc); % Primary Key (inline) increment (mode )
skey_inc_mode = mode(job_meta.skey_inc);
if pkey_inc_mode/pkey_inc_factor ~= floor(pkey_inc_mode/pkey_inc_factor) || ...
        skey_inc_mode/skey_inc_factor ~= floor(skey_inc_mode/skey_inc_factor)
    fprintf('Error: Interpolated grid has non-integer inline and crossline. \n')
else
    % Read seismic data
    [~, traces, ilxl_read, ~] = node_segy_read(job_meta_path,'1',i_block);
          
    pkey_inc_interp = pkey_inc_mode/pkey_inc_factor;
    skey_inc_interp = skey_inc_mode/skey_inc_factor;
    
%     block_keys_interp(:,1) = job_meta.block_keys(:,1)+job_meta.aperture;
%     block_keys_interp(:,2)  = job_meta.block_keys(:,2)-job_meta.aperture;
%     block_keys_interp(:,3)  = job_meta.block_keys(:,3)+job_meta.aperture;
%     block_keys_interp(:,4)  = job_meta.block_keys(:,4)-job_meta.aperture;

    % calculate new blocks keys
    block_keys_interp(:,1) = job_meta.block_keys(str2double(i_block),1)-pkey_inc_mode/2+job_meta.aperture;
    block_keys_interp(:,2) = job_meta.block_keys(str2double(i_block),2)+pkey_inc_mode/2-pkey_inc_interp-job_meta.aperture;
    block_keys_interp(:,3) = job_meta.block_keys(str2double(i_block),3)-skey_inc_mode/2+job_meta.aperture;
    block_keys_interp(:,4) = job_meta.block_keys(str2double(i_block),4)+skey_inc_mode/2-skey_inc_interp-job_meta.aperture;

%     block_keys_interp(:,1) = job_meta.block_keys(:,1)-pkey_inc_mode/2+job_meta.aperture;
%     block_keys_interp(:,2) = job_meta.block_keys(:,2)+pkey_inc_mode/2-pkey_inc_interp-job_meta.aperture;
%     block_keys_interp(:,3) = job_meta.block_keys(:,3)-skey_inc_mode/2+job_meta.aperture;
%     block_keys_interp(:,4) = job_meta.block_keys(:,4)+skey_inc_mode/2-skey_inc_interp-job_meta.aperture;

    pkey_min = min(ilxl_read(:,1));
    pkey_max = max(ilxl_read(:,1));
    skey_min = min(ilxl_read(:,2));
    skey_max = max(ilxl_read(:,2));
    skeyn = 1+((skey_max-skey_min)/skey_inc_mode);
    pkeyn = 1+((pkey_max-pkey_min)/pkey_inc_mode);
    
    % Create 3D volume
    n_iline = (ilxl_read(:,1)-pkey_min)/pkey_inc_mode+1;
    n_xline = (ilxl_read(:,2)-skey_min)/skey_inc_mode+1;
    lin_ind = ((n_iline-1).*skeyn)+n_xline;
    vol_3d = zeros(skeyn*pkeyn,job_meta.n_samples{1},'single');
    vol_3d(double(lin_ind),:) = traces';
    vol_3d = reshape(vol_3d,skeyn,pkeyn,job_meta.n_samples{1});
    vol_3d = permute(vol_3d,[3 1 2]);
    [x,y,z] = meshgrid(skey_min:skey_inc_mode:skey_max,0:1:job_meta.n_samples{1}-1,pkey_min:pkey_inc_mode:pkey_max);
    [xi,yi,zi] = meshgrid(skey_min:skey_inc_interp:skey_max,0:1:job_meta.n_samples{1}-1,pkey_min:pkey_inc_interp:pkey_max);
    
    clear ilxl_read
    il_i = squeeze(zi(1,:,:));
    xl_i = squeeze(xi(1,:,:));
    ilxl_read(:,1) = il_i(:); % new inlines
    ilxl_read(:,2) = xl_i(:); % new crosslines
    
    pkey_min = min(ilxl_read(:,1));
    pkey_max = max(ilxl_read(:,1));
    skey_min = min(ilxl_read(:,2));
    skey_max = max(ilxl_read(:,2));
    skeyn = 1+((skey_max-skey_min));
    pkeyn = 1+((pkey_max-pkey_min));
    
    ilxl_keep = ilxl_read(:,1) >= block_keys_interp(1) & ilxl_read(:,1) <= ...
        block_keys_interp(2) & ilxl_read(:,2) >= ...
        block_keys_interp(3) & ilxl_read(:,2) <= ...
        block_keys_interp(4);
    n_iline = (ilxl_read(ilxl_keep,1)-pkey_min)+1;
    n_xline = (ilxl_read(ilxl_keep,2)-skey_min)+1;
    lin_ind_keep = ((n_iline-1).*skeyn)+n_xline;
    
    ntraces = size(lin_ind_keep,1);
    resultno = 1;
    % Save outputs into correct structure to be written to SEGY.
    results_out{resultno,1} = 'Meta data for output files';
    results_out{resultno,2}{1,1} = ilxl_read(ilxl_keep,:);
    %results_out{resultno,2}{2,1} = uint32(zeros(size(traces{vol_index_wb},2),1));
    results_out{resultno,2}{2,1} = uint32(zeros(ntraces,1));
    
    ebdichdr = ['interp '];
    ebdichdr2{1,2} = '';
    tmpebc = '';
    ebcstrtowrite = sprintf('%-3200.3200s',[results_out{resultno,1} '  ' ebdichdr '  ' tmpebc]);
    results_out{resultno,1} = ebcstrtowrite;
    
    resultno = resultno + 1;
    results_out{resultno,1} = strcat(job_meta.volumes{1},'_interp');
    results_out{resultno,2} = interp3(single(x),single(y),single(z),vol_3d,single(xi),single(yi),single(zi));
    results_out{resultno,2} = reshape(results_out{resultno,2},job_meta.n_samples{1},[]);
    results_out{resultno,2} = results_out{resultno,2}(:,lin_ind_keep);
    results_out{resultno,3} = 0;
    resultno = resultno + 1;
    clear vol_3d
    
    % segy write function
    if exist(strcat(job_meta.output_dir,'interpolated_data/'),'dir') == 0
        output_dir = strcat(job_meta.output_dir,'interpolated_data/');
        mkdir(output_dir);
    else
        output_dir = strcat(job_meta.output_dir,'interpolated_data/');
    end
    
    i_block = str2double(i_block);
    node_segy_write(results_out,i_block,job_meta.s_rate/1000,output_dir);
end

end
