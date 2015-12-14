function calculate_3d_dip_st(job_meta_path,i_block,vol_index,start_slab,end_slab,sigma)

segy_plot_blocks(job_meta_path,'1');

scale_sigma = 4;
aperture = num2str(scale_sigma*str2num(sigma));

add_aperture_to_job_meta(job_meta_path,aperture)

job_meta = load(job_meta_path);  

pkey_inc_mode = mode(job_meta.pkey_inc); % Primary Key (inline) increment (mode )
skey_inc_mode = mode(job_meta.skey_inc);        

if str2double(end_slab) > job_meta.n_samples{str2double(vol_index)}
    end_slab = num2str(job_meta.n_samples{str2double(vol_index)});                    % update endslab if the data provided is of shorter length
end

n_slices = str2double(end_slab)-str2double(start_slab)+1;   % Number of Slices

%%
ebdichdr = ['dip parameters sigma ',num2str(sigma)];
ebdichdr2{1,2} = '';
tmpebc = '';

for ebcii = (size(ebdichdr2,1)-1):-1:1
    tmpebcc = regexp(ebdichdr2{ebcii,2},'/','split');
    tmpebc = [tmpebc tmpebcc{1}  tmpebcc{end}]; 
end
tmpebc = sprintf('%-3200.3200s',tmpebc);
clear tmpebcc ebdichdr2;

%%
loopfin = size(job_meta.liveblocks,1);                          % Number of live blocks
lpi = 348;

while lpi <= loopfin
    i_block = job_meta.liveblocks(lpi);                            % Block Number for Current Live Block
    [~, traces, ilxl_read, ~] = ...
        node_segy_read(job_meta_path,vol_index,num2str(i_block));
    traces = [zeros(1,size(traces,2)); traces(2:end,:)];
    
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
    vol_3d = zeros(skeyn*pkeyn,n_slices,'single');
    vol_3d(double(lin_ind),1:n_slices) = traces(str2double(start_slab):str2double(end_slab),:)';        
    
    vol_3d = reshape(vol_3d,skeyn,pkeyn,n_slices);
    vol_3d = permute(vol_3d,[3 1 2]);
    
    [dipx,dipy,planarity,cline,cplane,cfault] = calculate_st(vol_3d,sigma,scale_sigma);
    
    % Trim volumes to area without aperture
    il_min_keep = job_meta.block_keys(i_block,1)+str2num(aperture);
    il_max_keep = job_meta.block_keys(i_block,2)-str2num(aperture);
    xl_min_keep = job_meta.block_keys(i_block,3)+str2num(aperture);
    xl_max_keep = job_meta.block_keys(i_block,4)-str2num(aperture);
    
    ilxl_keep = ilxl_read(:,1) >= il_min_keep & ilxl_read(:,1) <= il_max_keep & ilxl_read(:,2) >= xl_min_keep & ilxl_read(:,2) <= xl_max_keep;
    n_iline = (ilxl_read(ilxl_keep,1)-pkey_min)/pkey_inc_mode+1;
    n_xline = (ilxl_read(ilxl_keep,2)-skey_min)/skey_inc_mode+1;
    lin_ind_keep = ((n_iline-1).*skeyn)+n_xline; 
    
    % Create attribute volumes 
    dipx = reshape(dipx,n_slices,skeyn*pkeyn);
    dipy = reshape(dipy,n_slices,skeyn*pkeyn);
    planarity = reshape(planarity,n_slices,skeyn*pkeyn);
    cline = reshape(cline,n_slices,skeyn*pkeyn);
    cplane = reshape(cplane,n_slices,skeyn*pkeyn);
    cfault = reshape(cfault,n_slices,skeyn*pkeyn);
    
    % Pad with zeros
    ntraces = size(lin_ind_keep,1);
    pad_zeros = zeros(str2double(start_slab)-1,ntraces);
    
    dipx = [pad_zeros; dipx(:,lin_ind_keep)];
    dipy = [pad_zeros; dipy(:,lin_ind_keep)];
    planarity = [pad_zeros; planarity(:,lin_ind_keep)];
    cline = [pad_zeros; cline(:,lin_ind_keep)];
    cplane = [pad_zeros; cplane(:,lin_ind_keep)];
    cfault = [pad_zeros; cfault(:,lin_ind_keep)];
    
    resultno = 1;
    % Save outputs into correct structure to be written to SEGY.
    results_out{resultno,1} = 'Meta data for output files';
    results_out{resultno,2}{1,1} = ilxl_read(ilxl_keep,:);
    %results_out{resultno,2}{2,1} = uint32(zeros(size(traces{vol_index_wb},2),1));
    results_out{resultno,2}{2,1} = uint32(zeros(ntraces,1));
    
    ebdichdr2{1,2} = '';
    tmpebc = '';
    
    ebcstrtowrite = sprintf('%-3200.3200s',[results_out{resultno,1} '  ' ebdichdr '  ' tmpebc]);
    results_out{resultno,1} = ebcstrtowrite;

    resultno = resultno + 1;
    results_out{1,3} = 'is_gather'; % 1 is yes, 0 is no

    results_out{resultno,1} = 'dipx';
    results_out{resultno,2} = dipx;
    results_out{resultno,3} = 0;
    resultno = resultno + 1;

    results_out{resultno,1} = 'dipy';
    results_out{resultno,2} = dipy;
    results_out{resultno,3} = 0;
    resultno = resultno + 1;

    results_out{resultno,1} = 'planarity';
    results_out{resultno,2} = planarity;
    results_out{resultno,3} = 0;
    resultno = resultno + 1;   
    
    results_out{resultno,1} = 'cline';
    results_out{resultno,2} = cline;
    results_out{resultno,3} = 0;
    resultno = resultno + 1; 
    
    results_out{resultno,1} = 'cplane';
    results_out{resultno,2} = cplane;
    results_out{resultno,3} = 0;
    resultno = resultno + 1; 
    
    results_out{resultno,1} = 'cfault';
    results_out{resultno,2} = cfault;
    results_out{resultno,3} = 0;
    resultno = resultno + 1; 
    
    node_segy_write(results_out,i_block,job_meta.s_rate/1000,[job_meta.output_dir,'/test/']);
    
    fprintf('-- Block %d of %d --\n',lpi,loopfin);
    lpi = lpi + 1; 
    clear  traces ilxl_read results_out
end

add_aperture_to_job_meta(job_meta_path,num2str(-1*str2num(aperture)));

end

function [dipx,dipy,planarity,cline,cplane,cfault] = calculate_st(vol_3d,sigma,scale_sigma)

sigma = str2num(sigma);

%% Calculate Structural Tensor
% Calculate gradients in x, y and z directions
[Ix, Iz, Iy] = gradient(vol_3d); % Ix crossline, Iz time, Iy inline 
% The first output FX is always the gradient
% along the 2nd dimension of F, going across columns.
% The second output FY is always the gradient along
% the 1st dimension of F, going across rows.  For
% the third output FZ and the outputs that follow,
% the Nth output is the gradient along the Nth dimension of F.

Ixy = Ix.*Iy;
Ixz = Ix.*Iz;
Iyz = Iy.*Iz;
Iyy = Iy.*Iy;
Ixx = Ix.*Ix;
Izz = Iz.*Iz;

clear Ix Iy Iz

% Smooth gradient components with a Gaussian
Ixy = imgaussian(Ixy,sigma,scale_sigma*sigma);
Ixz = imgaussian(Ixz,sigma,scale_sigma*sigma);
Iyz = imgaussian(Iyz,sigma,scale_sigma*sigma);
Iyy = imgaussian(Iyy,sigma,scale_sigma*sigma);
Ixx = imgaussian(Ixx,sigma,scale_sigma*sigma);
Izz = imgaussian(Izz,sigma,scale_sigma*sigma);

Ixy = Ixy(:);
Ixz = Ixz(:);
Iyz = Iyz(:);
Iyy = Iyy(:);
Ixx = Ixx(:);
Izz = Izz(:);

disp('Calculate slopes')
parfor i_loop = 1:length(Ixx);
    [l1,l2,l3,nx,ny,nz,~,~,~,~,~,~] = EigenVectors3D(Ixx(i_loop), Ixy(i_loop), Ixz(i_loop), Iyy(i_loop), Iyz(i_loop), Izz(i_loop));
    planarity(i_loop,1) = (l1-l2)./l1;
    cplane(i_loop,1) = (l1-l2)./(l1+l2);
    cline(i_loop,1) = (l2-l3)./(l2+l3);
    dipx(i_loop,1) = atand(nx/nz);
    dipy(i_loop,1) = atand(ny/nz);
end

cfault = (1-cplane).*cline;

end