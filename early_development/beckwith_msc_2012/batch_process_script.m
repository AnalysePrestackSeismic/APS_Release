% Prompt depending on algorithm chosen

% Start editing here %
input_dir = {'/data/TZA/dtect/bg_tza_site_survey_vol2_msc'}; % add as many directories as needed no slash at end
%/data/TZA/dtect/bg_tza_site_survey_vol2_msc

% paths to files instead?

% function call to scan all files in directories
files_in = directory_scan(input_dir);

% select files to input
index_files = input('Enter numbers of segy files to scan (in bracket [], e.g. [1 3 5]): ');
% this will used for all files
il_byte = 189;
xl_byte = 193;
extra_bytes_to_scan = [37];

% Do you want to scan all files that you have selected? [1 - Yes, 0 - No]
scan_files = 1;

files = length(index_files);

% scan the input files
input_segy = segy_read_files(scan_files,index_files,files_in,il_byte,xl_byte,extra_bytes_to_scan);

max_workers = 64; 
% Control the memory usage
n_blocks = 1;
output_dir = '/data/TZA/dtect/bg_tza_site_survey_vol2_msc/Misc/matlab_out';


n_freq = 110;
res_inc = 3;

% calculate number of traces per worker
input_segy{1}.aperture = 0;
%input_segy = segy_make_proc_pos(input_segy,aperture);

input_segy = njobs_worker(input_segy,max_workers,n_blocks,n_freq);
     
%matlabpool(max_workers)
sch = findResource('scheduler','configuration', defaultParallelConfig);
         
start_batch_index = 1;   

path_depend = {'/apps/gsc/matlab-mcode-beta/beckwith_msc_2012'};
    
for ii = 1:1:max_workers*n_blocks
     %fit_2d_gauss(filepath,n_samples,traces_process,n_freq,output_location)
     j{ii} = batch(sch,@fit_2d_gauss,0,{input_segy{1}.filepaths,input_segy{1}.n_samples,...
         input_segy{1}.trace_pointers(input_segy{1}.index_worker(ii,4):input_segy{1}.index_worker(ii,5),:),n_freq,res_inc,output_dir},'pathdependencies',path_depend,'matlabpool',0);            
      
     start_batch_index = start_batch_index + n_blocks;
end