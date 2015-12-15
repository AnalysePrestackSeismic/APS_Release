function [traces process_positions] = node_segy_read_traces(block_mat_all,block_mat,process_files_mat,index_file_read)
% need to add range of loading
    
    process_files = load(process_files_mat,'-mat');
    process_positions = node_make_processing_positions(block_mat_all,block_mat);
    
    fid = open_segy_file(cell2mat(strcat(process_files.path(index_file_read),process_files.name(index_file_read))));

    if strcmp(process_positions.distribute_type,'slice') && isfield(process_files, 'index_wb')    
        wb_diff_index = max(process_files.index_wb(:,3)) - min(process_files.index_wb(:,3));
        n_samples = wb_diff_index+(max(process_positions.z_pos)-min(process_positions.z_pos)+1);
        if n_samples > process_positions.n_samples
            n_samples = process_positions.n_samples;
        end
    else
        n_samples = max(process_positions.z_pos)-min(process_positions.z_pos)+1;
    end
    
    % we always read in regular manner so let's work what we have read
    [~,I_min] = min(process_positions.ilxl_pos(:,3));   
    [~,I_max] = max(process_positions.ilxl_pos(:,3));  

    [~, Locb_min] = ...
        ismember(process_positions.ilxl_pos(I_min,1:2),process_positions.trace_ilxl_bytes(:,1:2), 'rows'); 
    [~, Locb_max] = ...
        ismember(process_positions.ilxl_pos(I_max,1:2),process_positions.trace_ilxl_bytes(:,1:2), 'rows');

    read_traces = process_positions.trace_ilxl_bytes(Locb_min:Locb_max,1:2);
    n_traces = Locb_max-Locb_min+1; 
    
    fseek(fid,process_positions.start_byte,'bof'); 
    skip = ((process_positions.n_samples-n_samples)*4)+process_positions.trc_head_length; % +(process_positions.n_samples*process_positions.ilxl_step)*4;
    prec = strcat(num2str(n_samples),'*uint32=>uint32');
    temp = fread(fid,[n_samples n_traces],prec,skip);
       
    if isfield(process_files, 'index_wb')
            
            % compare what we have read with what we want
            [Lia2, Locb2] = ...
            ismember(process_positions.ilxl_pos(:,1:2),read_traces(:,1:2), 'rows');
            
            if size(temp,2) > size(process_positions.ilxl_pos,1)
                temp = ibm2double(temp(:,nonzeros(Locb2))); 
                logical_size = 1;
                % index wb compare
                [~, Locb3] = ...
                ismember(read_traces(:,1:2),process_files.index_wb(:,1:2), 'rows');
                index_wb = process_files.index_wb(nonzeros(Locb3),3);     
            else
                temp = ibm2double(temp);
                logical_size = 0;
                [~, Locb3] = ...
                ismember(read_traces(:,1:2),process_files.index_wb(:,1:2), 'rows');
                index_wb = process_files.index_wb(nonzeros(Locb3),3);       
            end                
            
            if strcmp(process_positions.distribute_type,'trace')
                % index wb compare         
                shift_temp_traces = NaN(size(temp,1),size(temp,2));   
                constant = 1; % just for display at the moment should be 1 and then trim matrix j    
                for ii = 1:1:size(temp,2)
                    shift_temp_traces(1:size(temp,1)-index_wb(ii)+constant,ii) = ...
                        temp(index_wb(ii):end,ii);
                end 
                shift_n_samples = n_samples - min(process_files.index_wb(:,3));
                [~,nan_index] = max(process_positions.z_pos == shift_n_samples);
            else
                shift_temp_traces = NaN(size(temp,1),size(temp,2));   
                constant = 1; % just for display at the moment should be 1 and then trim matrix j    
                index_wb = index_wb - min(process_files.index_wb(:,3)) + 1;
                
                for ii = 1:1:size(temp,2)
                    shift_temp_traces(1:size(temp,1)-index_wb(ii)+constant,ii) = ...
                        temp(index_wb(ii):end,ii);
                end
                nan_index = length(process_positions.z_pos); 
            end
            clearvars temp
            traces = NaN(nan_index,length(process_positions.ilxl_pos));
            
            z = process_positions.z_pos-min(process_positions.z_pos)+1;
            z = z(1:nan_index);
            process_positions.z_pos = z;
            if logical_size == 1
                traces(:,Lia2) = shift_temp_traces(z,:);
            else
                traces(:,Lia2) = shift_temp_traces(z,nonzeros(Locb2));   
            end
    else % no wb flatten
            % let's compare what we have read to what we want
            [Lia2, Locb2] = ...
            ismember(process_positions.ilxl_pos(:,1:2),read_traces(:,1:2), 'rows');
            % maximum size of array to fill
            traces = NaN(length(process_positions.z_pos),length(process_positions.ilxl_pos));
            z = process_positions.z_pos-min(process_positions.z_pos)+1;            
            traces(:,Lia2) = ibm2double(temp(z,nonzeros(Locb2)));
    end

%     positions = NaN(length(traces),2);
%     positions(Lia2,1:2) = read_traces(nonzeros(Locb2),1:2);   

    close_segy_file(fid);    
    
end

function data2=ibm2double(data)
% Convert IBM 32-bit floating-point format to IEEE 64-bit floating-point format
%      (based on the algorithm of function "ibm2num" written by Brian Farrelly)
% See also: ibm2single
%
% Written by: E. Rietsch: July 19, 2009 
% Last updated:
%
%         data2=ibm2double(data)
% INPUT 
% data    a matrix of unsigned 4-byte integers (uint32) representing data
%         in IBM floating-point format
% OUTPUT
% data2   corresponding matrix of doubles (IEEE format)

% get sign from first bit
% get exponent from first byte, last 7 bits
% remove bias from exponent 
% get mantissa from last 3 bytes

data2=(1-2*double(bitget(data,32))).*16.^ ...
  (double(bitshift(bitand(data,uint32(hex2dec('7f000000'))),-24))-64) .* ...
  (double(bitand(data,uint32(hex2dec('00ffffff'))))/2^24);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fid = open_segy_file(filepaths)
   
    mformat='ieee-be';
    fid=fopen(filepaths,'r',mformat);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function close_segy_file(fid)
   
    fclose(fid);

end