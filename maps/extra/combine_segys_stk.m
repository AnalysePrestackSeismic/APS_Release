function [] = combine_segys_stk(job_meta_path,inputpath)
% -------------------------------------------------------------------------
% merge all the seperate stack volumes
% -------------------------------------------------------------------------


job_meta = load(job_meta_path);             % Load job meta information 

if job_meta.is_gather == 0                  % For the case of using angle stacks
    i_vol_max = job_meta.nvols;             % Number of angle stacks
else
    i_vol_max = 0;
end

%==================================================================================================
% now scan the input directory and gat all the files to merge

[files_in,nfiles] = directory_scan(filepath,filename_string); % files_in is a structure, .names cell array, .path cell array 
                                                              % directory_scan filters out the file names with the user supplied string and also returns the number of such files: nfiles that ...  
                                                              % ... can be used for index pre alocation in future steps 
files_in.names = sort_nat(files_in.names);                    % Natural Sort file names 
start_point = pwd;                                            % Remember starting directory

% need to sort the files names into block order to write the data into the
% output segy file
for i_file = 1:1:nfiles
        
        filename = files_in.names{i_file};                          % Sequentially put file name in the structure
        filepath = files_in.path{i_file}; 

end

%

    % check segy write functions - many different versions now!
    if exist(strcat(job_meta.output_dir,'bg_combined_angle_stks/'),'dir') == 0
        output_dir = strcat(job_meta.output_dir,'bg_combined_angle_stks/');
        mkdir(output_dir);
    else
        output_dir = strcat(job_meta.output_dir,'bg_combined_angle_stks/');
    end

loopfin = size(job_meta.liveblocks,1);      % Number of Live Blocks
lpi = 2;                                    % Loop Index
count = 1;                                  % Counter for blocks that have a corresponding wavelet file

%====================================================================================================
% first loop to open the output file and write the ebcdic and binary header

% read first dataset in
i_block = job_meta.liveblocks(1);     % Reference the block numbers for live blocks
tmpfileout = sprintf('opening file %s',strcat(inputpath,'fft_wavelets_block_',num2str(i_block),'.bin'));
  
% read all the data for this block
% node_segy_read(job_meta_path,vol_index,i_block)
[~, vol_traces, ilxl_read, offset_read] = node_segy_read(job_meta_path,'1',i_block);










n_results = size(results_out,1);
for i_results = 2:1:n_results % top row in results is meta inforation includes output directory (can either be local to node or on system)
    
    
    
    
    
    
    
    %file_name = fullfile(output_dir, sprintf('%s_result_block_%d%s',results_out{i_results,1},i_block,'.segy'));
    %file_name = fullfile(output_dir, sprintf('%s_block_%d',results_out{i_results,1},i_block),'.segy');
    out_dir = strcat(output_dir,results_out{i_results,1},'/');
    
    if ~exist(out_dir,'dir')
       mkdir(out_dir); 
    end
    
    file_name = fullfile(out_dir, sprintf('%s_block_%d%s',results_out{i_results,1},i_block,'.segy'));
    
    % separate folder for each i_results
    
    % Define Variables
    % start_time = 0;
    % sample_rate = 4;
    n_samples = size(results_out{i_results,2},1);
    
%     if n_samples < n_samples_out
%        pad_samples = n_samples_out - n_samples;
%        results_out{i_results,2} = [results_out{i_results,2}; zeros(size(results_out{i_results,2},2),pad_samples)];        
%        n_samples = n_samples_out;
%     end
    
    coscal_val = -100; % Coordinate scalar to be written to SEGY.
    block_sz = 100; % Number of traces written at a time
    
    %% Initialise files and headers
    fid_ilxl_f32 = fopen(file_name,'w'); % Open file for writing
    header = zeros(60,block_sz,'int32'); % Define zero matrix in dimensions of header.
    
    % Set an array for the coordinate scalar as int16 packed into int32 data
    tmpcolscal = zeros(1,(block_sz*2),'int16');
    tmpcolscal(1:2:size(tmpcolscal,2)) = coscal_val;
    tmpcolscal2 = typecast(tmpcolscal,'int32');
    
    % Define binary header
    bin_header = zeros(200,1,'uint16'); % Define zero vector as size of trace header
    bin_header(7) = 1;                  % Number of data traces per ensemble
    bin_header(9) = sample_rate*1000;   % Sample rate (needs to be in microseconds)
    bin_header(11) = n_samples;         % Number of samples per trace
    bin_header(13) = 5;                 % SEGY data number format 5 = ieee float 754
    bin_header(14) = 1;                 % Ensemble fold
    bin_header(15) = 4;                 % Trace sorting code 4 = horizontally stacked , 2 = cdp ensemble
    bin_header(28) = 1;                 % Measurement system
    bin_header(151) = 256;              % SEGY revision number
    bin_header(152) = 1;                % Fixed length flag  1= size always the same
    
    % Write out the ebcdic_header
    % seismic.text_header = ['AAAAA',blanks(3190),'ZZZZZ'];
    seismic.text_header = sprintf('%-3200.3200s',results_out{1,1});
    fwrite(fid_ilxl_f32,seismic.text_header,'char*1',0,'ieee-be');
    
    % Write out the binary header
    fwrite(fid_ilxl_f32,bin_header,'uint16',0,'ieee-be');
    












while lpi <= loopfin
    i_block = job_meta.liveblocks(lpi);     % Reference the block numbers for live blocks
    tmpfileout = sprintf('opening file %s',strcat(job_meta.wav_directory,'fft_wavelets_block_',num2str(i_block),'.bin'));
  
% read all the data for this block
% node_segy_read(job_meta_path,vol_index,i_block)
[~, vol_traces, ilxl_read, offset_read] = node_segy_read(job_meta_path,'1',i_block);
    
    
    
    


    resultno = 1;
    % Save outputs into correct structure to be written to SEGY.
    results_out{resultno,1} = 'Meta data for output files';
    results_out{resultno,2}{1,1} = ilxl_read;
    results_out{resultno,3} = 'is_gather'; % 1 is yes, 0 is no
    results_out{resultno,2}{2,1} = uint32(zeros(size(angle_stk{aidx},2),1));
    %was written as uint32(zeros(ntraces,1));
    %results_out{resultno,2}{2,1} = offset_read';
     
    ebcstrtowrite = sprintf('%-3200.3200s',[results_out{resultno,1} '  ' ebdichdr '  ' tmpebc]);
    results_out{resultno,1} = ebcstrtowrite;
    
    resultno = resultno + 1;
    
    % correct file names added by SRW - 02/07/14
    
%     if kk == startvol
%         testdiscpt = ['gathers_segy_io_',num2str(origstart_vol),'_',num2str(volinc),'_',num2str(origstart_vol+angwidth)];
%     else
        testdiscpt = ['angle_stk_range_',num2str(kk),'_',num2str((kk+(angwidth-volinc)))];
%     end
       
    results_out{resultno,1} = strcat(testdiscpt);
    %results_out{2,2} = digi_intercept;
    results_out{resultno,2} = angle_stk{aidx};
    results_out{resultno,3} = 0;
    aidx = aidx +1;
    

    
    
    
end
    
    
end