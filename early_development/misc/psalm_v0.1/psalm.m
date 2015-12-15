% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% PSALM - Parallel Seismic Analysis Leveraging MATLAB
% Date: 03 October 2012
% Authors: James Selvage and Jonathan Edgar
%
%%%%%% NEED AUTOMATIC PATH UPDATING FOR ADDING NEW FUNCTIONS
% Function to run
[distribute_type algo_name func_name] = algorithm_to_run();
process_files.distribute_type = distribute_type;
process_files.algo_name = algo_name;
process_files.func_name = func_name;

input_dir = {input('Enter input data directory: ')};
process_files.path_for_blocks = input('Enter output directory for node job files: ');

[files_in_dir] = directory_scan(input_dir);
index_files = ...
    input('Enter numbers of segy files to scan (in bracket [], e.g. [1 3 5]): ');
process_files.nfiles = length(index_files);

% Scan through each file and ask about type
for i_file = 1:1:process_files.nfiles
    process_files.path(i_file) = files_in_dir.path(index_files(i_file));
    process_files.name(i_file) = files_in_dir.names(index_files(i_file));
    if i_file == 1
        % Display the type of segy files that the algorithms can handle
        segy_types();
    end
    process_files.type_code(i_file) = ...
        input(sprintf('Enter file type for %s: ',process_files.name{i_file}));  
    process_files.type_description(i_file) = ...
        segy_types(process_files.type_code(i_file));    
    if strcmp(process_files.type_description{i_file},'Angle stack');
        process_files.angle(i_file) = ...
            input(sprintf('  File is angle stack, enter angle (degrees): '));
    else 
        process_files.angle(i_file) = NaN;
    end
end

% We can flatten on the water bottom on the fly; parting the ocean will be
% in version 2.
process_files.in_flatten_water_bottom = input('\nDo you want to flatten on the water bottom? (1 - Yes, 0 - No): ');
%in_flatten_water_bottom = input;

if process_files.in_flatten_water_bottom == 1
    if sum(process_files.type_code == 1) >= 1 % then we have a full stack...halle
        [~,match] = max((process_files.type_code == 1));
        index_to_scan = match;
        process_files.index_wb_scan = index_to_scan;
    elseif sum(process_files.type_code == 1) == 0 && ...
            sum(process_files.type_code == 2) >= 1 % then we have angle stacks
        index_to_scan = process_files.angle == ...
            min(process_files.angle(process_files.angle > 20));
        process_files.index_wb_scan = index_to_scan;
    elseif sum(process_files.type_code == 1) == 0 ...
            && sum(process_files.type_code == 2) == 0 
        fprintf('\nI could try flattening using this volume,\n but conversion cannot always be performed.\n');
        in_flatten_water_bottom = 0;
        index_to_scan = 1;        
    end
else
    index_to_scan = 1;
    process_files.index_wb_scan = index_to_scan;
end

% We need to make a structure for one of the files - who will be the chosen
% one? Either it should be the first file, or a full stack or middle angle 
% stack if flattening on water bottom is desired
bytes = input(sprintf('Enter inline and crossline byte locations for %s [inline crossline]: ',...
    process_files.name{index_to_scan}));
process_files.ilbyte = bytes(1);
process_files.xlbyte = bytes(2);
process_files.extra_bytes = input('Enter extra bytes to scan [byte_1 byte_2]: ');
seismic = segy_make_structure(process_files.path{index_to_scan},...
   process_files.name{index_to_scan},process_files.ilbyte,...
   process_files.xlbyte,process_files.extra_bytes);
process_files.binary_header = seismic.binary_header;
process_files.n_samples = seismic.n_samples;
process_files.s_rate = seismic.s_rate;
process_files.n_iline = seismic.n_iline;
process_files.n_xline = seismic.n_xline; 
process_files.il_inc = seismic.il_inc;
process_files.xl_inc = seismic.xl_inc;
process_files.min_iline = seismic.min_iline;
process_files.max_iline = seismic.max_iline;
process_files.min_xline = seismic.min_xline;
process_files.max_xline = seismic.max_xline;

% Save the process_files structure so that it is available to the nodes
process_files_mat = strcat(process_files.path_for_blocks,func_name,'_process_files.mat');
save(process_files_mat,'-struct','process_files','-v7.3');

if process_files.in_flatten_water_bottom == 1
    wb_track.name = process_files.name(index_to_scan);
    wb_track.path = process_files.path(index_to_scan);
    wb_track.ilbyte = bytes(1);
    wb_track.xlbyte = bytes(2);
    wb_track.n_blocks = 1000;
    wb_track.path_for_blocks = process_files.path_for_blocks;
    wb_track.func_name = 'water_bottom_flatten';
    wb_track.func_to_run = process_files.func_name;
    wb_track.n_samples = seismic.n_samples;
    wb_track.s_rate = seismic.s_rate;
    wb_track.n_iline = seismic.n_iline;
    wb_track.n_xline = seismic.n_xline; 
    wb_track.min_iline = seismic.min_iline;
    wb_track.max_iline = seismic.max_iline;
    wb_track.min_xline = seismic.min_xline;
    wb_track.max_xline = seismic.max_xline;
    wb_track.il_inc = seismic.il_inc;
    wb_track.xl_inc = seismic.xl_inc;
    wb_track.processing_grid = segy_make_processing_grid(1,0,0,1,0,wb_track);

    wb_mat = strcat(wb_track.path_for_blocks,wb_track.func_name,'_process_files.mat');
    save(wb_mat,'-struct','wb_track','-v7.3');
    batch_distribute_processing_grid(seismic,wb_track);

    fprintf('\n*** Tracking the water bottom ... ***\n\n');

%     Compile function
%     compile_function(wb_track.func_name);
%     
%     % System call to SLURM
%     fprintf('Submitting task to SLURM');
%     % if varargin
%     node_slurm_submit(wb_track.func_name,wb_track.n_blocks,wb_track.path_for_blocks,wb_track.n_blocks);
    
    % need to add whiles to functions for them to wait for SLURM processes to
    % finish
    
%     for i_block = 1:1:wb_track.n_blocks
%         block_mat_all = strcat(wb_track.path_for_blocks,'water_bottom_flatten','_positions_block_all.mat');
%         block_mat = strcat(wb_track.path_for_blocks,'water_bottom_flatten','_positions_block_',num2str(i_block),'.mat');
% 
%         % water_bottom_flatten('/home/thedailyrant/node_jobs/test3/water_bottom_flatten_positions_block_1.mat','/home/thedailyrant/node_jobs/test3/water_bottom_flatten_process_files.mat')
%         water_bottom_flatten(block_mat_all,block_mat,wb_mat,num2str(wb_track.n_blocks))
%     end 
end

% check structures are the same in other files

% make the processing grid
fprintf('File %s grid is:\n',process_files.name{index_to_scan}); 
fprintf('Inline step: %d\n',seismic.il_inc); 
fprintf('Crossline step: %d\n',seismic.xl_inc); 
fprintf('Z sample rate (ms): %d\n',seismic.s_rate/1000); 
fprintf('Z samples: %d\n',seismic.n_samples); 

fprintf('\nEnter the processing grid that you want to use:\n');
ilxl_step = input('Enter inline/crossline step: ');
ilxl_aperture = input('Enter inline/crossline aperture: ');
ilxl_aperture_step = input('Enter inline/crossline aperture step: ');
z_step = input('Enter z step: ');
z_aperture = input('Enter z aperture: ');

process_files = load(process_files_mat,'-mat');

% specify amount to parallelise
process_files.n_blocks = input('Enter number of blocks to divide processing into: ');
process_files.processing_grid = segy_make_processing_grid(ilxl_step,ilxl_aperture,...
    ilxl_aperture_step,z_step,z_aperture,process_files);
save(process_files_mat,'-struct','process_files','-v7.3');

batch_distribute_processing_grid(seismic,process_files);

% fprintf('\n*** Compiling the function ***\n\n');
% 
% compile_function(process_files.func_name);
% 
% fprintf('\n\n*** Distributing the compiled function using SLURM ***\n\n');
% %anomalous_body_connector(block_mat_all,block_mat,process_files_mat,anomalous_threshold,connectivity,join_block_id)
% node_slurm_submit(process_files.func_name,process_files.n_blocks,process_files.path_for_blocks,0.995,6,5);
% %node_slurm_submit(process_files.func_name,process_files.n_blocks,process_files.path_for_blocks,process_files.path_for_blocks);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %