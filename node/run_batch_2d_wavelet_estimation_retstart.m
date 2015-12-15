function[]=run_batch_2d_wavelet_estimation_retstart(algo_run_qc_path,algorithm_name,slurm_part,n_cores,decimate,first_live_block)

batch_run_file_info=load(algo_run_qc_path);

for g=1:batch_run_file_info.nfiles
    id=strfind(batch_run_file_info.algo_run_qc{g},'wckey=');
     batch_run_file_info.wkey_orig{g}=batch_run_file_info.algo_run_qc{g}((id+6):end);
 
    batch_run_file_info.algo_run_qc_restart{g}=node_slurm_submit2013_restart_CAN(algorithm_name,batch_run_file_info.job_meta_files{g},slurm_part,n_cores,batch_run_file_info.wkey_orig{g},decimate,first_live_block);
    fprintf('submitted line : %g \n',g)
    pause(5);
    
    
end

save(algo_run_qc_path,'-struct','batch_run_file_info','-v7.3'); % Saves Seismic structure to mat file


end