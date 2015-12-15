%function [] = sas_test_script(seismic_mat_path,i_block,n_blocks,window_length,binary_slice_path,block_path,wb_path,output_dir)
algorithm_name = 'seismic_anomaly_spotter_lite_sw';
binary_slice_path = '/data/Global/dtect/SRW_test_dataset/SEGY/matlab/full_run/digi_results/';
n_blocks = 139;
seismic_mat_path = '/data/Global/dtect/SRW_test_dataset/SEGY/matlab/full_run/20-25_angle_stack_SRW.mat_lite';


all_windows = [1; 5; 10; 50; 100; 150; 200; 250; 500];
output_directory{1,1} = '/data/Global/dtect/SRW_test_dataset/SEGY/matlab/full_run/sas_results/window_1/';
output_directory{2,1} = '/data/Global/dtect/SRW_test_dataset/SEGY/matlab/full_run/sas_results/window_5/';
output_directory{3,1} = '/data/Global/dtect/SRW_test_dataset/SEGY/matlab/full_run/sas_results/window_10/';
output_directory{4,1} = '/data/Global/dtect/SRW_test_dataset/SEGY/matlab/full_run/sas_results/window_50/';
output_directory{5,1} = '/data/Global/dtect/SRW_test_dataset/SEGY/matlab/full_run/sas_results/window_100/';
output_directory{5,1} = '/data/Global/dtect/SRW_test_dataset/SEGY/matlab/full_run/sas_results/window_150/';
output_directory{5,1} = '/data/Global/dtect/SRW_test_dataset/SEGY/matlab/full_run/sas_results/window_200/';
output_directory{5,1} = '/data/Global/dtect/SRW_test_dataset/SEGY/matlab/full_run/sas_results/window_250/';
output_directory{5,1} = '/data/Global/dtect/SRW_test_dataset/SEGY/matlab/full_run/sas_results/window_500/';


for test = 1:length(all_windows)
    window_length = num2str(all_windows(test))
    %output_dir = output_directory{test,1}
    output_dir = output_directory{1,1}
    
    %for i_block = 1:n_blocks;
        node_slurm_submit_lite_sw(algorithm_name,n_blocks,seismic_mat_path,window_length,binary_slice_path,output_dir);
        %seismic_anomaly_spotter_lite_sw(seismic_mat_path,i_block,n_blocks,window_length,binary_slice_path,output_dir);
    %end
    
    
%     sas_to_segy_sw(seismic_mat_path,i_block,n_blocks,block_path,wb_path,output_dir)
    
    
end 
