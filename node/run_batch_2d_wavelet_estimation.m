function[]=run_batch_2d_wavelet_estimation(input_filepath,filename_string,output_filepath,il_byte,xl_byte,offset_byte,anggath,algorithm_name,slurm_part,n_cores,decimate,first_live_block)

run_scan_input=1;
parallel=1;
%% ----------------scan for files----------------
[files_in,nfiles] = directory_scan(input_filepath,filename_string); % files_in is a structure, .names cell array, .path cell array 
                                                              % directory_scan filters out the file names with the user supplied string and also returns the number of such files: nfiles that ...  
                                                              % ... can be used for index pre alocation in future steps 
%   files_in.names = sort_nat(files_in.names);                    % Natural Sort file names 
%start_point = pwd;                                            % Remember starting directory
% batch_run_file_info = load('/data/CAN/segy/NE_Newfoundland_Flemish_Pass/BG_internal_production/digi/wavelet_estimationandQC_CAN_run_info.mat');
batch_run_file_info.files_in=files_in;
batch_run_file_info.nfiles=nfiles;
algo_run_qc_path='/data/CAN/segy/NE_Newfoundland_Flemish_Pass/BG_internal_production/digi_2/wavelet_estimationandQC_CAN_run_info.mat';

%% -------------make folders-------------------------------------------------------------------

if exist(output_filepath,'dir') == 0
    mkdir(output_filepath)
    system(['chmod 777 ',output_filepath]);
end

links_dir=strcat(output_filepath,'links/');

if exist(links_dir,'dir') == 0
    mkdir(links_dir)
    system(['chmod 777 ',links_dir]);
end

results_dir=strcat(output_filepath,'result/');


if exist(results_dir,'dir') == 0
    mkdir(results_dir)
    system(['chmod 777 ',results_dir]);
end

%% --------------------make links and output directories----------------------------------------------------

for f=1:nfiles
    
    file_name_temp=strcat(files_in.path{f},'/',files_in.names{f});
    id=strfind(files_in.names{f},'.');
    links_input_dir{f}=strcat(links_dir,files_in.names{f}(1:(id-1)),'/');
    batch_run_file_info.links_input_dir{f}=links_input_dir{f};
    
    if exist(links_input_dir{f},'dir') == 0
        mkdir(links_input_dir{f})
        system(['chmod 777 ',links_input_dir{f}]);
    end
    link_name_temp=strcat(links_input_dir{f},'gather_',algorithm_name,'_input_link_block1.segy');
    system(['ln -s ',file_name_temp,' ',link_name_temp]);
    
    output_dir{f}=strcat(results_dir,files_in.names{f}(1:(id-1)),'/');
    
    if exist(output_dir{f},'dir') == 0
        mkdir(output_dir{f})
        system(['chmod 777 ',output_dir{f}]);
    end
    
    
    
end

clear f file_name_temp id  link_name_temp ;

    


%% ------------------- run segy_make job----------------------

if(run_scan_input==1)
    for g=1:nfiles
        job_meta_files{g} = segy_make_job(links_input_dir{g},'block',il_byte,xl_byte,offset_byte,parallel,anggath,output_dir{g});
        fprintf('Completed segy make job for %d out of %d lines \n',g,nfiles);
        batch_run_file_info.job_meta_files{g}=job_meta_files{g};
        batch_run_file_info.algo_run_qc{g}=node_slurm_submit2013_CAN(algorithm_name,job_meta_files{g},slurm_part,n_cores,decimate,first_live_block);
        %     pause(10);
        fprintf('submitted line : %g \n',g)
    end
end
%% ----------------------run algo--------------------------------------------

%algo_run_qc_path=strcat(output_filepath,algorithm_name,'_run_info.mat');
% batch_run_file_info = struct;
%batch_run_file_info = load('/data/CAN/segy/NE_Newfoundland_Flemish_Pass/BG_internal_production/FK_Gathers_mica/trim_calculation_can_run_info.mat');
% 
% batch_run_file_info.links_input_dir=links_input_dir;
% batch_run_file_info.job_meta_files=job_meta_files;

% for g=15:nfiles
%     
%     
%     
% end

save(algo_run_qc_path,'-struct','batch_run_file_info','-v7.3'); % Saves Seismic structure to mat file


end