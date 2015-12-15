function [] = sw_segy_write_clean_tg(results_out,i_block, n_blocks, sample_rate, output_dir, varargin)
%% Function to write SEGY files.
% The "results_out" format is a cell array which always contains metadata in the
% first row. Subsequent rows contain data matricies to be written. This
% function loops through, writing each to a seperate SEGY file.


n_results = size(results_out,1);
for i_results = 2:1:n_results % top row in results is meta inforation includes output directory (can either be local to node or on system)
    %file_name = fullfile(output_dir, sprintf('%s_result_block_%d%s',results_out{i_results,1},i_block,'.segy'));
    %file_name = fullfile(output_dir, sprintf('%s_block_%d',results_out{i_results,1},i_block),'.segy');
    file_name = fullfile(output_dir, sprintf('%s_block_%d%s',results_out{i_results,1},i_block,'.segy'));
    %% Define Variables
    % start_time = 0;
    % sample_rate = 4;
    n_samples = size(results_out{i_results,2},1);
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
    bin_header(152) = 1;                % Fixed length flag  1 = size always the same
    
    % Write out the ebcdic_header
    seismic.text_header = ['AAAAA',blanks(3190),'ZZZZZ'];
    fwrite(fid_ilxl_f32,seismic.text_header,'char*1',0,'ieee-be');
    
    % Write out the binary header
    fwrite(fid_ilxl_f32,bin_header,'uint16',0,'ieee-be');
    
    %% Loop through the data and write out to SEGY
    tic
    for ii= 1:block_sz:size(results_out{i_results,2},2)
        
        maxblock = ii + block_sz - 1;
        
        if maxblock >  size(results_out{i_results,2},2)
            maxblock = size(results_out{i_results,2},2);
            block_sz = (maxblock - ii+1);
            header = zeros(60,block_sz,'int32');
            tmpcolscal = zeros(1,(block_sz*2),'int16');
            tmpcolscal(1:2:size(tmpcolscal,2)) = coscal_val;
            tmpcolscal2 = typecast(tmpcolscal,'int32');
        end
        
        temparr = typecast(single(reshape(results_out{i_results,2}(:,ii:maxblock),1,(n_samples*(maxblock-ii+1)))),'int32');
        temparr2 = [header; reshape(temparr,n_samples,block_sz)];
        
        % Set header values
        count = ii:1:maxblock;
        temparr2(1,:) = count;
        
        % Set byte 29 to 1 as a 16 bit integer into a int32, so the first 16
        % bytes have to represent 1 when read as unit16, so multiply 1 by 2^16
        % which is 65536
        temparr2(8,:) = 65536;
        
        % To write the int16 coordinate scalar as -100 as int32 in to byte location 71
        temparr2(18,:) = tmpcolscal2;
        
        % write the inline crossline numbers to bytes 189 and 193
        temparr2(48,:) = typecast(results_out{1,2}{1,1}(ii:maxblock,1)','int32'); % assumes inline / crossline are the same for all files
        temparr2(49,:) = typecast(results_out{1,2}{1,1}(ii:maxblock,2)','int32');
        
        fwrite(fid_ilxl_f32,temparr2,'int32',0,'ieee-be');
    end
    fclose(fid_ilxl_f32);
    
    toc
    
end