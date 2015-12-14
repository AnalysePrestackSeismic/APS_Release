function [] = node_segy_write(results_out,i_block, sample_rate, output_dir, varargin)
%% ------------------ Disclaimer  ------------------
% 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) makes no representation or warranty, express or implied, in 
% respect to the quality, accuracy or usefulness of this repository. The code
% is this repository is supplied with the explicit understanding and 
% agreement of recipient that any action taken or expenditure made by 
% recipient based on its examination, evaluation, interpretation or use is 
% at its own risk and responsibility.
% 
% No representation or warranty, express or implied, is or will be made in 
% relation to the accuracy or completeness of the information in this 
% repository and no responsibility or liability is or will be accepted by 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) in relation to it.
%% ------------------ License  ------------------ 
% GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
%% github
% https://github.com/AnalysePrestackSeismic/
%% ------------------ FUNCTION DEFINITION ---------------------------------
% node_segy_write: function to read traces from a specific block with a
% scanned segy volume
%   Arguments:
%       results_out = a cell array which always contains metadata in the
%           first row. Subsequent rows contain data matrices to be written
%       i_block = integer specifying which block results are from
%       sample_rate = sample of SEGY to add to segy header
%       output_dir = path specifying output location
%       varargin = does not do anything at the moment
%
%   Outputs:
%       none
%
%   Writes to Disk:
%       Loops through entries in results_out and creates separate SEGY
%       files at location specified by output_dir

%%
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
    
    %% Define Variables
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
    
    %% Loop through the data and write out to SEGY
    last_byte = 0;
    
    
    bytes_per_sample = 4;
    skip_textual_binary = 3600;
    trc_head = 240;
    trc_length = n_samples*bytes_per_sample;
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
        
        % Set header values to byte loc 1 as 4 byte integer
        count = ii:1:maxblock;
        temparr2(1,:) = count;
        
        % Set byte 29 to 1 as a 16 bit integer into a int32, so the first 16
        % bytes have to represent 1 when read as unit16, so multiply 1 by 2^16
        % which is 65536
        temparr2(8,:) = 65536;
        
        % To write the int16 coordinate scalar as -100 as int32 in to byte location 71
        temparr2(18,:) = tmpcolscal2;
        
        % to write the no of samples into trace header 115 as 16 bit integer
        temparr2(29,:) = n_samples;
        
        % to write the sample rate into trace header 117 as 16 bit integer
        temparr2(30,:) = sample_rate*1000*65536;   % Sample rate (needs to be in microseconds)
        
        % write the inline crossline numbers to bytes 189 and 193
        %temparr2(48,:) = cast(results_out{1,2}{1,i_results-1}(ii:maxblock,1)','int32'); % assumes inline / crossline are the same for all files
        %temparr2(49,:) = cast(results_out{1,2}{1,i_results-1}(ii:maxblock,2)','int32');
        %temparr2(48,:) = typecast(results_out{1,2}{1,i_results-1}(ii:maxblock,1)','int32'); % assumes inline / crossline are the same for all files
        %temparr2(49,:) = typecast(results_out{1,2}{1,i_results-1}(ii:maxblock,2)','int32');
        
        % the test of n_results is to cover the ability to write results
        % sets with different sets of trace headers (pkey skey tkey) they
        % would go as more rows into resultsout{1,2} and would the looped
        % through, rather gave up on this and tend to call multiple calls
        % to node_segy_write
        
        if size(results_out{1,2},2) < n_results
            % Write offset information to byte 37 as 4 byte 32 bit integer
            temparr2(10,:) = typecast(results_out{1,2}{2,1}(ii:maxblock,1)','int32');
            % write the inline crossline numbers to bytes 189 and 193
            temparr2(48,:) = results_out{1,2}{1,1}(ii:maxblock,1)';
            temparr2(49,:) = results_out{1,2}{1,1}(ii:maxblock,2)';
        else
            % Write offset information to byte 37 as 4 byte 32 bit integer
            temparr2(10,:) = typecast(results_out{1,2}{2,i_results-1}(ii:maxblock,1)','int32');
            % write the inline crossline numbers to bytes 189 and 193
            temparr2(48,:) = results_out{1,2}{1,i_results-1}(ii:maxblock,1)';
            temparr2(49,:) = results_out{1,2}{1,i_results-1}(ii:maxblock,2)';
        end
        
        fwrite(fid_ilxl_f32,temparr2,'int32',0,'ieee-be');
        
        % make byte location array % cj this needs to be a 64 bit int
        byte_locs(ii:maxblock,1) = last_byte+(trc_head:trc_head+trc_length:block_sz*(trc_length+trc_head));
        last_byte = byte_locs(end,1)+trc_length; % store to add on during loop                                    
    end
    fclose(fid_ilxl_f32);
    
    % file has written successfully so make .lite file for future use
    % should be able to just use gather version of compress
    byte_locs(:,1) = byte_locs(:,1)+skip_textual_binary;  

    % see above for discussion on this if
    if size(results_out{1,2},2) < n_results
        if results_out{i_results,3} == 0
            compress_ilxl_bytes = trace_compress_ilxl_bytes([results_out{1,2}{1,1} byte_locs],size(results_out{i_results,2},2));
        else
            compress_ilxl_bytes = gather_compress_ilxl_bytes([results_out{1,2}{1,1} byte_locs results_out{1,2}{2,1}],size(results_out{i_results,2},2));
        end        
    else
        if results_out{i_result,3} == 0
            compress_ilxl_bytes = trace_compress_ilxl_bytes([results_out{1,2}{1,1} byte_locs],size(results_out{i_results,2},2));
        else
            compress_ilxl_bytes = gather_compress_ilxl_bytes([results_out{1,2}{1,i_results-1} byte_locs results_out{1,2}{2,i_results-1}],size(results_out{i_results,2},2));
        end
    end
    clear byte_locs;
%             if anggath == 1
%                 compress_ilxl_bytes = gather_compress_ilxl_bytes(trace_ilxl_bytes,blocktr);
%             else
%                 compress_ilxl_bytes = gather_compress_ilxl_bytes_offset(trace_ilxl_bytes,blocktr);
%             end
    
    
    
    file_mat = fullfile(out_dir, sprintf('%s_block_%d%s',results_out{i_results,1},i_block,'.mat_orig_lite'));
    
    filepath_binary = uint64(file_name);
    pad_filepath = zeros(1,(2000-length(filepath_binary)));
    filepath_binary = [filepath_binary,pad_filepath];
    fid_write = fopen(file_mat,'w');        
    %fwrite(fid_write,[filepath_binary';seismic.file_type;seismic.s_rate;seismic.n_samples;seismic.n_traces;il_byte;xl_byte;offset_byte],'double');  
    % need to set is_gather
    fwrite(fid_write,filepath_binary','double'); 
    fwrite(fid_write,[bin_header(13);bin_header(9);bin_header(11);size(results_out{i_results,2},2);189;193;37;results_out{i_results,3}],'double'); 
    fwrite(fid_write,reshape(compress_ilxl_bytes',[],1),'double');  
    fclose(fid_write);
    
end

end