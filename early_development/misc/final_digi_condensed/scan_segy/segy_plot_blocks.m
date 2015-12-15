function segy_plot_blocks(job_meta_path,i_vol)

job_meta = load(job_meta_path);

[seismic, ~, ~, ~] = node_segy_read(job_meta_path,i_vol,'1');

figure
scatter(seismic.trace_ilxl_bytes(:,1),seismic.trace_ilxl_bytes(:,2))
hold all
scatter(seismic.trace_ilxl_bytes(:,1),seismic.trace_ilxl_bytes(:,4))

    for i_block = 1:1:size(job_meta.block_keys,1)
        colour = [random('bino',1,0.5),random('bino',1,0.5),random('bino',1,0.5)];
        fill([job_meta.block_keys(i_block,2);job_meta.block_keys(i_block,2);...
            job_meta.block_keys(i_block,1);job_meta.block_keys(i_block,1)],...
            [job_meta.block_keys(i_block,3);job_meta.block_keys(i_block,4);...
            job_meta.block_keys(i_block,4);job_meta.block_keys(i_block,3)],colour)
    end

end