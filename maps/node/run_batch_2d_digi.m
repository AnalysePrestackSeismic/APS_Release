function[]=run_batch_2d_digi(algo_run_qc_path,algorithm_name,slurm_part,n_cores,startvol,volinc,endvol,tottracerun,maxzout,wavevar)

%% Program to run DIGI on a series of 2D lines in a folder

batch_run_file_info=load(algo_run_qc_path); % Load the batch run infor mat file

%   Loop through all the lines
for g=1:batch_run_file_info.nfiles
    
    wavelet_avg(batch_run_file_info.job_meta_files{g});         % Run wavelet Averaging
    meta_data_2d_smo_xy(batch_run_file_info.job_meta_files{g}); % Run smoothening of the meta data
    batch_run_file_info.digi_run_qc{g}=node_slurm_submit2013_CAN(algorithm_name,batch_run_file_info.job_meta_files{g},slurm_part,n_cores,startvol,volinc,endvol,tottracerun,maxzout,wavevar);   % Run DIGI
    fprintf('submitted line : %g \n',g)
    pause(5);
    
    
end

save(algo_run_qc_path,'-struct','batch_run_file_info','-v7.3'); % Save the batch run infor mat file
close all;

%print all the sacct commands to see job running sttus
for k=1:52
    fprintf(batch_run_file_info.digi_run_qc{k});
end

end