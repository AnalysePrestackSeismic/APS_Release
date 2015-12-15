%% Example chain

% Run psalm
psalm();

% Wavelet estimation
for qq = 1:6
    block_mat = sprintf('/usr/local/data/node_jobs/tza_test/wavelets/wavelet_estimation_positions_block_%d.mat',qq);
    process_files_mat = '/usr/local/data/node_jobs/tza_test/wavelets/wavelet_estimation_process_files.mat';
    output_dir = '/usr/local/data/mat/tza_test/wavelets/';
    wavelet_estimation(block_mat, process_files_mat, output_dir);
    %save(sprintf('chi_%d',qq),'chi','-v7.3');
end

% Minimum energy chi
for qq = 1:6
    block_mat = sprintf('/usr/local/data/node_jobs/tza_test/chi_model/minimum_energy_chi_positions_block_%d.mat',qq);
    process_files_mat = '/usr/local/data/node_jobs/tza_test/chi_model/minimum_energy_chi_process_files.mat';
    output_dir = '/usr/local/data/mat/tza_test/chi_model/';
    minimum_energy_chi(block_mat, process_files_mat, output_dir);
    %save(sprintf('chi_%d',qq),'chi','-v7.3');
end

% Invert AVA
for qq = 89:250
    block_mat=sprintf('/usr/local/data/node_jobs/tza_test2/int_grad_inv_proj_positions_block_%d.mat',qq);
    process_files_mat='/usr/local/data/node_jobs/tza_test2/int_grad_inv_proj_process_files.mat';
    wavelet_mat_dir='/usr/local/data/mat/tza_test/wavelets/';
    chi_mat_dir='/usr/local/data/mat/tza_test/chi_model/';
    output_dir='/usr/local/data/mat/tza_test2/';
    
    invert_for_ava_v6(block_mat, process_files_mat, wavelet_mat_dir, chi_mat_dir, output_dir);
    
    fprintf('Completed block %d of %d\n',qq,250)
end

% Save segy
segy_write_traces(block_mat);

% Run psalm again
psalm();

% SAS test
for qq = 100:106
    block_mat=sprintf('/usr/local/data/mat/tza_test/sas/seismic_anomaly_spotter_positions_block_%d.mat',qq);
    process_files_mat='/usr/local/data/node_jobs/tza_test/sas/seismic_anomaly_spotter_process_files.mat';
    output_dir='/usr/local/data/mat/tza_test/sas/';
     
    seismic_anomaly_spotter(block_mat, process_files_mat,output_dir);
end

% Save segy
segy_write_traces(block_mat);

% Anomaly connect
block_path = '~/node_jobs/test10/';
process_path = '~/node_jobs/test10/anomalous_body_connector_process_files.mat';

for i_block = 1:1:n_block;
    
    block_cat = strcat(block_path,'anomalous_body_connector_positions_block_',num2str(i_block),'.mat');
    anomalous_body_connector(block_cat,process_path,0.95,6,n_block); 
    
end

% Anomaly merge
merge_anomalous_body_connector('anomalous_body_connector','~/node_jobs/test_connect10/')

% Anomaly report
report_anomalous_body_connector(block_path,10);

% Save segy
segy_write_traces(block_path);


