function[]=check_run_batch_2d_(algo_run_qc_path)

batch_run_file_info=load(algo_run_qc_path);
counter=1;

for g=1:batch_run_file_info.nfiles
    
    job_meta=load(batch_run_file_info.job_meta_files{g});
    
    if length(job_meta.liveblocks)<24
        
    batch_run_file_info.lower_liveblocks{counter}=batch_run_file_info.files_in.names{g};
    fprintf(batch_run_file_info.files_in.names{g});fprintf('\n');
    counter=counter+1;
    end
end

for g=1:batch_run_file_info.nfiles
    
    batch_run_file_info.algo_run_qc{g}
    fprintf('\n');
end

save(algo_run_qc_path,'-struct','batch_run_file_info','-v7.3');
    
end