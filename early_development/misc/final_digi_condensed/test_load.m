function test_load(i_block)

job_meta_path = '/data/NAM/segy/Project_Oryx/3d_stacks/output/job_meta/job_meta_29Nov2013.mat';
job_meta = load(job_meta_path);



[~, traces, ilxl_read] = node_segy_read(job_meta_path,'1',i_block);

results_out{1,1} = 'Meta data for output files';
results_out{1,2}{1,1} = ilxl_read;
results_out{1,2}{2,1} = uint32(zeros(size(traces,2)));
results_out{2,1} = 'test_load_write';
results_out{2,2} = traces;

node_segy_write_traces2(results_out,str2double(i_block),5,job_meta.output_dir)
end