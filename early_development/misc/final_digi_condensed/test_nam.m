%%
%segy_make_job({'/data/NAM/dtect/pro_oy/Misc/'},'block','189','193','37','0','/data/NAM/dtect/pro_oy/Misc/','10','1','0')
%%
jobmetapath = '/data/NAM/segy/Project_Oryx/3d_stacks/job_meta/job_meta_12Nov2013.mat';
i_vol = '1';
for i_block = 1:1:101 
    job_meta = load(jobmetapath);

    [seismic, traces, ilxl_read] = node_segy_read(jobmetapath,i_vol,num2str(i_block));
    [wb_idx] = water_bottom_flatten_lite(traces);
    [traces] = trace_flatten(traces,wb_idx);
    
    flat_directory = strcat(job_meta.output_dir,'slices/');
    mkdir(flat_directory);
    % Write slice ordered binary of the minimum energy projection for use in SAS
    fid_out = fopen(strcat(flat_directory,'flat_full_stack','_slices_block_',num2str(i_block),'.bin'),'w');
    nsamples = size(traces,1);
    fwrite(fid_out,nsamples,'float32');
    fwrite(fid_out,traces','float32');
    fclose(fid_out);
    
end
%%
% for i_slab = 1:1:200
%     seismic_anomaly_spotter(jobmetapath,'slices/','1',num2str(i_slab),'200','100')
%     %seismic_anomaly_spotter(job_meta_path,slice_file_rel_path,i_vol,i_slab,n_slab,window_length)
% end