function [] = digi
% DIGI: Top level program to run Dynamic Intercept Gradient Inversion, the
% Seismic Anomaly Spotter, and Anomalous Body Connector programs along with
% all ancillary functions. Guides the user through running the various
% functions on their dataset using command line prompts.
fprintf('Welcome to PSALM lite: \n\n Step 1: SEGY geometry scan: \n\n Please follow the instructions to scan the geometry of your SEGY file: \n');

% Collect inputs:
% filepath = '/data/CHN/segy/6316/BG6316/ANGLESTACKFINAL/DIGI/input_anglestk_from_dtect/'; % input('Please enter input data directory containing the angle stacks you wish to use(in single quotes): ');
% file_to_scan = '35-45.sgy';                 % input('Please enter the filename of the SEGY file you wish to scan(in single quotes): ');
% il_byte = 189;                                          % input('Please enter the inline byte location: ');
% xl_byte = 193;                                          % input('Please enter the crossline byte location: ');
% output_dir = '/data/scratch/wrigley/'; %input('Please enter the directory in which you would like all result and temporary files to be saved: ');

filepath = '/data/Global/dtect/SRW_test_dataset/SEGY/'; % input('Please enter input data directory containing the angle stacks you wish to use(in single quotes): ');
filename = '20-25_angle_stack_SRW.sgy';                 % input('Please enter the filename of the SEGY file you wish to scan(in single quotes): ');
il_byte = 189;                                          % input('Please enter the inline byte location: ');
xl_byte = 193;                                          % input('Please enter the crossline byte location: ');
output_dir = '/data/Global/dtect/SRW_test_dataset/SEGY/matlab/digi_condensed/'; %input('Please enter the directory in which you would like all result and temporary files to be saved: ');

% filepath = input('\nPlease enter input data directory containing the angle stacks you wish to use(in single quotes, followed by a forward slash / ): ');
% file_to_scan = input('\nPlease enter the filename of the SEGY file you wish to scan(in single quotes): ');
% il_byte = input('\nPlease enter the inline byte location: ');
% xl_byte = input('\nPlease enter the crossline byte location: ');
% output_dir = input('\nPlease enter the output directory in which you would like all result files to be saved (in single quotes, followed by a forward slash / ): ');

% filepath = '/segy/NOR/VBT1_pl599_2013/2012_preSDM_deliverables/final_deliverables/final_angle_volumes/'%input('\nPlease enter input data directory containing the angle stacks you wish to use(in single quotes, followed by a forward slash / ): ');
% file_to_scan = '20_25DEG_SUB_ANGLE_STACK_TIME.segy'%input('\nPlease enter the filename of the SEGY file you wish to scan(in single quotes): ');
% il_byte = 189%input('\nPlease enter the inline byte location: ');
% xl_byte = 193 %input('\nPlease enter the crossline byte location: ');
% output_dir =  '/segy/NOR/VBT1_pl599_2013/2012_preSDM_deliverables/final_deliverables/IG_deliverables_SRW/with_wavelet_normalisation/' %input('\nPlease enter the output directory in which you would like all result files to be saved (in single quotes, followed by a forward slash / ): ');

% filepath = '/data/CHN/segy/cgh_DIGI/SAS_Inputs/'; % input('Please enter input data directory containing the angle stacks you wish to use(in single quotes): ');
% file_to_scan = 'EER_Me_crop.segy';                 % input('Please enter the filename of the SEGY file you wish to scan(in single quotes): ');
% il_byte = 189;                                          % input('Please enter the inline byte location: ');
% xl_byte = 193;                                          % input('Please enter the crossline byte location: ');
% output_dir = '/data/CHN/segy/cgh_DIGI/SAS_Outputs/'; %input('Please enter the directory in which you would like all result and temporary files to be saved: ');


% Create log file containing processing flow.
log_fid = fopen(strcat(output_dir,'psalm_log.txt'),'a');
fprintf(log_fid,'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
% fprintf(log_fid,strcat('PSALM run initiated: ', date,'\n\nInput Directory:',filepath,'\nOutput Directory:',output_dir,'\n\nAngle Stacks Used:\n'));
fprintf(log_fid,sprintf('PSALM_lite run initiated: %s \n\nInput Directory: %s\nOutput Directory: %s\n\nAngle Stacks Used:\n',datestr(clock),filepath,output_dir));


% Scan directory for files
start_point = pwd;
cumm_files = 1;
cd(filepath);
[~,nfiles] = (system('ls -B *gy | wc -l'));
nfiles=str2double(nfiles);

[~,fnames] = system('ls -B1 *gy');
numeric = double(fnames);

% Preallocate memory for some variables
count = 2;
fname_index = zeros(1,nfiles+1);

% Loop to separate out each file name from one long character string
for ij= 1:length(fnames)
    if numeric(1,ij) == 10
        fname_index(1,count) = ij;
        count = count+1;
    end
end

% Loop to read each file as ascii into a cell array
for ik=1:nfiles
    files_in.names{cumm_files+ik-1} = fnames(1,(fname_index(1,ik)+1):(fname_index(1,ik+1)-1));
    % files_in.path{cumm_files+ik-1} = filepath{ii};
    fprintf(log_fid,strcat(files_in.names{cumm_files+ik-1},'\n'));
end

fclose(log_fid);

fprintf('\nThe following files have been found:\n')
for il = 1:nfiles
    angle_range{il} = str2double(regexp(files_in.names{il},'(\d{2})','match'));
    fprintf('File %d: %s, Angle range: %d-%d \n',il,files_in.names{il},angle_range{il});
end

mid_stk = floor(nfiles/2);

corr_angles = input('\nAre the angle ranges listed above correct? Please enter y or n: \n','s');
if strcmpi(corr_angles,'n');
    file_to_change = input('Please enter the file numbers which require changing in brackets [] separated by commas e.g. [1,3,4]: \n')
    
    for aa = 1:length(file_to_change)
        bb = file_to_change(aa);
        angle_range{bb} = input(sprintf('Please enter the correct angle range for file %d in brackets [] separated by a comma e.g. [10,15]: ',bb));
    end
    
    fprintf(['\n Thank you. Proceeding to scan SEGY geometry...\n'])
else
    fprintf(['\n Proceeding to scan SEGY geometry...\n'])
end

%file_to_scan = files_in.names{mid_stk}; % Finds the mid stack to get geometry & water bottom from.
fprintf(sprintf('Now scanning mid angle stack: %s to obtain SEGY geometry and waterbottom...\n',filename))

[seismic] = segy_make_structure(filepath, filename, il_byte, xl_byte, output_dir);


for kk = 1:1:nfiles
    angle_stacks{kk,1} = files_in.names{kk};        % Picks up filename
    angle_stacks{kk,2} = angle_range{kk};
    angle_stacks{kk,3} = angle_stacks{kk,2}(1) + ((angle_stacks{kk,2}(2) - angle_stacks{kk,2}(1))/2); % finds the mid angle
end
seismic.angles = [angle_stacks{:,3}];
seismic.angle_range = angle_range;
seismic.filepath = filepath;
seismic.input_filenames = {angle_stacks(:,1)};
seismic.nfiles = nfiles;

% Save extra info to mat file
file_mat = dir(strcat(seismic.output_dir,'segy_structure/','*.mat'));
file_mat = regexp(file_mat.name,'\.','split');
save(strcat(seismic.output_dir,'segy_structure/',file_mat{1}),'-struct','seismic','-v7.3');
seismic_mat_path = strcat(seismic.output_dir,'segy_structure/',file_mat{1},'.mat');


fprintf(['\nSEGY scan complete. \n \n Step 2: Wavelet estimation \n Please '...
    'follow the onscreen prompts to estimate time varying wavelets'...
    ' from your input SEGY file \n\n'])

n_blocks = num2str((input('Please enter the number of blocks to divide the volume into for processing:')));
i_block = '1';

%write info to log file
log_fid = fopen(strcat(output_dir,'psalm_log.txt'),'a');
fprintf(log_fid,strcat('\nNumber of blocks used for computation: ',n_blocks,'\n'));
fclose(log_fid);


% Check for existing wavelets and generate if none found
if ((exist(strcat(seismic.output_dir,'wavelets/'),'dir')) && (length(dir(strcat(seismic.output_dir,'wavelets/','*.mat'))) == (seismic.nfiles*2)+2))
    w = input(['\n It appears wavelets have already been estimated for these angle stacks. \n Please press:'...
        '\n \n 1 - to continue with the existing wavelets \n 2 - to re-calculate the wavelets\n']);
    if w == 2
        
        wavelet_avg(seismic_mat_path,n_blocks);
        
    end
else
    fprintf(sprintf('\n Now calculating time varying wavelets: %d jobs for each %d angle stacks \n',str2num(n_blocks),nfiles));
    %n_blocks = str2num(n_blocks);
    wavelet_avg(seismic_mat_path,n_blocks);
end




% Section to run DIGI and/or SAS
algorithm = input(['\n Step 3: Please enter the number of the algorithm you wish '...
    'to run: \n \n 1 = Dynamic Intercept/Gradient Inversion'...
    '\n 2 = Seismic Anomaly Spotter\n\n\']);
if algorithm == 1
    nodes_to_use = input('\n Preparing to run DIGI...\n\n In the CPU monitor window which appears, please observe which nodes are free prior to running any jobs on the cluster');
    
    system('/apps/gsc/scripts/clustercpu_hal all');
    
    n_blocks = str2double(n_blocks);
    node_slurm_submit_lite_sw('int_grad_inv_proj_lite_sw',n_blocks,seismic_mat_path)
    
    
else if algorithm == 2
        if ((exist(strcat(seismic.output_dir,'digi_results/','minimum_energy_slices/'),'dir') ~= 0) && (length(dir(strcat(seismic.output_dir,'digi_results/','minimum_energy_slices/','digi_minimum_energy_eer_projection_slices_block_*.bin'))) == str2double(n_blocks)))
            fprintf('\n %d minimum energy slices found. Submitting %d slices to SAS for processing...\n',str2double(n_blocks), str2double(n_blocks))
            window_length = num2str(input('\nPlease enter the window length you would like to use (50 is default): '));
            if exist(strcat(seismic.output_dir,'sas_results'),'dir') == 0
                mkdir(strcat(seismic.output_dir,'sas_results'));
            end
            n_blocks = str2double(n_blocks);
            node_slurm_submit_lite_sw('seismic_anomaly_spotter_lite_sw',n_blocks,seismic_mat_path,window_length)
            
        else
            aa = input(['\nNo minumum energy slices have been found. Please press 1 to return'...
                'and run DIGI or 2 if you wish to run SAS on a full stack volume: \n']);
            if aa == 1
                fprintf('/n Preparing to run DIGI\n');
                node_slurm_submit_lite_sw('int_grad_inv_proj_lite_sw',n_blocks,seismic_mat_path)
            else if aa ==2
                    fprintf('\n This section of code has not been added yet...')
                end
            end
        end
    end
end



%% Alternative section to compile the entire job then run it.
% algorithm = input(['\nIn order to submit your job to the cluster, please enter which algorithms you wish to run\n'...
%     'in square brackets seperated by commas in the order in which they are to be submitted. \nFor example; to run DIGI followed by SAS, and then \n'...
%     'the anomalous body connector enter [1,2,3] \n\n Algorithms available: \n\n1 - Dynamic Intercept Gradient Inversion\n'...
%     '2 - Seismic Anomaly Spotter\n3 - Anomalous Body Connector\n\n'])
%
% for zz = 1:1:size(algorithm)
%     if algorithm(zz) == 1
%         algorithm_name = 'int_grad_inv_proj_lite_sw'
%         total_expected_files = n_blocks*3;
%     else if algorithm(zz) == 2
%             algorithm_name = 'seismic_anomaly_spotter_lite_sw'
%             total_expected_files =
%         else if algorithm(zz) == 3
%                 fprintf('Algorithm not yet implimented')
%             end
%         end
%     end
%
%     node_slurm_submit_lite_sw(algorithm_name,n_blocks,seismic_mat_path)
%
%     wait_flag = 0;
%     while wait_flag <= total_expected_files
%         [~,wait_flag] = system('ls -B * | wc -l');
%         wait_flag = str2double(wait_flag);
%
%         if wait_flag == total_expected_files;
%             wait_flag = total_expected_files +1;
%         end
%     end
%
% end




% add loop to submit to cluster and then count files till ready to
% progress to next task.l


end
