% scale velocity volume from Vvert to Vnmo using delta volume

% process volume in inlines so change block size in meta data

% - run segy_make_job for velocity and delta volumes
% - find max number of xls
% - change block size so that each block consists of complete inlines
% - read each block from each volume
% - Vnmo = Vvert x sqrt(1 + 2 x delta)

%---------------------------------------------------------------------------------
%% scan velocity volume

vel_pathname = '/data/TZA/segy/2012_kussini_pgs_full_sequence/presdm/velocities/';

vel_outdir = '/data/TZA/segy/2012_kussini_pgs_full_sequence/presdm/velocities/vel_scan/'

vel_file = 'Kussini_Vint_depth_VTI';
vnmo_file = 'Kussini_Vint_depth_Isotropic';

segy_make_job(strcat(vel_pathname,'velocity_link/'),vel_file,'189','193','37','0','0',vel_outdir);

[vel_meta_files,nfiles] = directory_scan({strcat(vel_outdir,'/job_meta/')},'job_meta_');

vel_meta_file = strcat(vel_meta_files.path{1},vel_meta_files.names{1});

vel_meta = load(vel_meta_file);

%---------------------------------------------------------------------------------
%% scan delta voiume

delta_pathname = '/data/TZA/segy/2012_kussini_pgs_full_sequence/presdm/velocities/';

delta_outdir = '/data/TZA/segy/2012_kussini_pgs_full_sequence/presdm/velocities/delta_scan/';

delta_file = 'Kussini_Delta_depth_VTI';

segy_make_job(strcat(vel_pathname,'delta_link/'),delta_file,'189','193','37','0','0',delta_outdir);

[delta_meta_files,nfiles] = directory_scan({strcat(delta_outdir,'/job_meta/')},'job_meta_');

delta_meta_file = strcat(delta_meta_files.path{1},delta_meta_files.names{1});

delta_meta = load(delta_meta_file);

%---------------------------------------------------------------------------------
%% change block size so each block is number of complete inlines
%

num_vel_xls = (1 + (vel_meta.skey_max - vel_meta.skey_min) / vel_meta.skey_inc);

vel_inline_size = num_vel_xls * vel_meta.n_samples{1} * 4;

max_file_size = 1024^3;

num_inlines_per_block = floor(max_file_size/vel_inline_size);

segy_block_size(vel_meta_file,0,num_inlines_per_block,num_vel_xls);

segy_block_size(delta_meta_file,0,num_inlines_per_block,num_vel_xls); % Assume delta file is same dimensions as velocities

vel_meta = load(vel_meta_file);
delta_meta = load(delta_meta_file);


%---------------------------------------------------------------------------------
%% loop through the blocks doing the vel conversion
%

for block = 1:vel_meta.n_blocks
    [vel_hdr,vel_traces,vel_ilxl,vel_offsets] = node_segy_read(vel_meta_file,'1',num2str(block));
    [delta_hdr,delta_traces,delta_ilxl,delta_offsets] = node_segy_read(delta_meta_file,'1',num2str(block));
    
    if ~isequal(delta_ilxl,vel_ilxl)
        disp('Inline/crossline ranges from delta and velocity volumes not identical');
        exit
    end
    
    vnmo = vvert2vnmo(vel_traces,delta_traces);
    
    segy_out{1,1} = 'Meta data for output files';
    segy_out{1,2}{1,1} = vel_ilxl;
    segy_out{1,2}{2,1} = int32(zeros(vel_hdr.n_traces,1));
    segy_out{1,1} = sprintf('%-3200.3200s','Vnmo calculated from Vvert and delta');
    segy_out{1,3} = 'is_gather'; % 1 is yes, 0 is no
    segy_out{2,1} = vnmo_file;
    segy_out{2,2} = vnmo;
    segy_out{2,3} = 0;
    
    node_segy_write(segy_out,block, vel_hdr.s_rate/1000, vel_pathname); % remember to divide srate by 1000 when writing back
    
end

segy_make_job(strcat(vel_pathname,vnmo_file,'/'),vnmo_file,'189','193','37','0','0',strcat(vel_pathname,vnmo_file,'/'));




%%---------------------------------------------------------------------------------
    
    