function [] = seismic_anomaly_spotter_traces(job_meta_path,vol_index,i_slab,n_slab,window_length)
% -------------------------------------------------------------------------
% SEISMIC_ANOMALY_SPOTTER: Takes the minimum energy projection from the
% dynamic intercept gradient inversion. Forms a cumulative distribution of
% the data which is used to write out a volume of data probability which
% can then be thresholded and fed into the anomalous body connector.
%   Inputs:
%       seismic_mat_path = path of metadata .mat file.
%       i_block = current block to be processed.
%       n_blocks = total number of blocks to submit for processing.
%       window_length = dictates the nuber of binary slices to be used for
%       CDF calculation
%       
% -------------------------------------------------------------------------

% Load job meta information 
job_meta = load(job_meta_path);

start_slice_idx = (1+((str2double(i_slab)-1)*floor(job_meta.n_samples{str2double(vol_index)}...
    /str2double(n_slab))))-str2double(window_length);

% Calculate the last trace to read for this slab (adjusts for extra slices in last slab)
if str2double(i_slab) == str2double(n_slab)
    end_slice_idx = window_length+job_meta.n_samples{str2double(vol_index)};
else
    end_slice_idx = str2double(window_length)+(str2double(i_slab)*floor(job_meta.n_samples{str2double(vol_index)}/str2double(n_slab)));
end

pad_start = [];
pad_end = [];
% vol_traces needs to be padded size
% Correct for window overshoots
if start_slice_idx <= 0
    pad_start = zeros(job_meta.vol_traces(str2double(vol_index)),1+abs(start_slice_idx));
    %pad_start = zeros(seismic.n_traces/n_blocks,1+abs(start_slice_idx));
    start_slice_idx = 1;
end
if end_slice_idx >= job_meta.n_samples{str2double(vol_index)}
    % if end_slice_idx >= ns
    pad_end = zeros(job_meta.vol_traces,end_slice_idx-job_meta.n_samples{str2double(vol_index)});
    end_slice_idx = job_meta.n_samples{str2double(vol_index)};
    %     end_slice_idx = ns;
end

% Calculate the number of slices to read in this slab
n_slices_to_read = end_slice_idx-start_slice_idx+1;
% Pre-allocate memory for the output zero-padded data probability slab
pkeyn = 1+((job_meta.pkey_max-job_meta.pkey_min)/job_meta.pkey_inc);
skeyn = 1+((job_meta.skey_max-job_meta.skey_min)/job_meta.skey_inc);
% generate volume for storing slices
slices = zeros(job_meta.vol_traces(str2double(vol_index)),n_slices_to_read);
% Read traces for each block
% Extract time samples and populate slices slab
for i_block = 1:1:job_meta.n_blocks
    [~, traces, ilxl_read] = ...
        node_segy_read(job_meta_path,num2str(vol_index),num2str(i_block));
    
    % Loop round to pick the water bottom
    [wb_idx{i_block}] = water_bottom_picker(traces,0); 
    
    % Position these as linear indices
    n_iline = ilxl_read(:,1)-job_meta.pkey_min+1;
    n_xline = ilxl_read(:,2)-job_meta.skey_min+1;
    lin_ind = ((n_iline-1).*skeyn)+n_xline;
    
    % Trim traces
    traces = traces(start_slice_idx:end_slice_idx,:)';
    
    slices(double(lin_ind),:) = traces;
    
end
slices = [pad_start,cell2mat(slices),pad_end];
n_windows = size(slices,2)-2*str2double(window_length);
dataprob = zeros(pkeyn*skeyn,n_windows);


% Loop through compressed file format to calculate linear indices
% Load seismic compress to get geometry
vol_index = str2double(vol_index);
i_row_count = 1;
for i_files = 1:1:size(job_meta.files,1)
    if strfind(job_meta.files{i_files},job_meta.volumes{vol_index})
        seismic = segy_read_binary(strcat(job_meta.paths{1},job_meta.files{i_files}));
        seismic.trace_ilxl_bytes = seismic.trace_ilxl_bytes(:,[1:2,4:end]);  
        % Calculate linear indices
        constant = 10000;        
        pkeys = unique(seismic.trace_ilxl_bytes(:,1));
        
        for i_row = 1:1:size(pkeys,1)                                 
            il_rows = find(seismic.trace_ilxl_bytes(:,1) == pkeys(i_row,1)); 
            
            if size(il_rows,1) > 1
                sub_trace_ilxl_bytes = seismic.trace_ilxl_bytes(il_rows,:);
                
                % check that crossline are consistent of padding is
                % required                
                for ii_xl = 1:1:size(il_rows,1)-1
                    if sub_trace_ilxl_bytes(ii_xl+1,2) - sub_trace_ilxl_bytes(ii_xl,3) == sub_trace_ilxl_bytes(ii_xl,4)
                        % no padding is required                          
                        sub_trace_ilxl_bytes(ii_xl,3) = sub_trace_ilxl_bytes(ii_xl+1,3);
                        sub_trace_ilxl_bytes(ii_xl+1,2) = sub_trace_ilxl_bytes(ii_xl,2);                        
                        % remove updated row
                    end                    
                end
                trace_ilxl_bytes(i_row_count,:) = unique(sub_trace_ilxl_bytes,'rows');
            else
                trace_ilxl_bytes(i_row_count,:) = seismic.trace_ilxl_bytes(il_rows,:);
            end   
            i_row_count = i_row_count + 1;
        end
        i_row_count = i_row_count + 1;
    end
end
 
for i_row = 1:1:size(trace_ilxl_bytes,1)
    cdplbl_min = trace_ilxl_bytes(i_row,1)*(job_meta.skey_min+constant);
    cdplbl_max = trace_ilxl_bytes(i_row,1)*(job_meta.skey_max+constant); 
    
    cdplbl_start = trace_ilxl_bytes(i_row,1)*(trace_ilxl_bytes(i_row,2)+constant); 
    if cdplbl_start > cdplbl_min
        n_traces_to_pad_before(i_row,1) = (cdplbl_start-cdplbl_min)/trace_ilxl_bytes(i_row,1);
    else
        n_traces_to_pad_before(i_row,1) = 0;
    end
    cdplbl_end = trace_ilxl_bytes(i_row,1)*(trace_ilxl_bytes(i_row,3)+constant);
    if cdplbl_end < cdplbl_max
        % need to pad
        n_traces_to_pad_after(i_row,1) = (cdplbl_max-cdplbl_end)/trace_ilxl_bytes(i_row,1);
    else
        n_traces_to_pad_after(i_row,1) = 0;
    end
    n_traces_row(i_row,1) = (cdplbl_end - cdplbl_start)/trace_ilxl_bytes(i_row,1);
end
n_traces_row = n_traces_row+1;
n_traces_to_pad_before = n_traces_to_pad_before+1;
live_traces = zeros(size(trace_ilxl_bytes,1),2);
live_traces(1,:) = [n_traces_to_pad_before(1),n_traces_row(1)+n_traces_to_pad_before(1)-1];
for i_row = 2:1:size(trace_ilxl_bytes,1)    
   live_traces(i_row,:) = [n_traces_to_pad_after(i_row-1)+n_traces_to_pad_before(i_row);n_traces_row(i_row)+n_traces_to_pad_before(i_row)+n_traces_to_pad_after(i_row-1)-1];
end



%Read the slices for this slab
[slices] = node_binary_read_slices(job_meta_path,slice_file_rel_path,num2str(vol_index),i_slab,n_slab,window_length);
%imagesc(reshape(slices(:,200),xln,iln))
% Work out the number of windows to loop through
n_windows = size(slices,2)-2*str2double(window_length);

dataprob = zeros(pkeyn*skeyn,n_windows);
window_length = str2double(window_length);
for ii = 1:n_windows
    % Calculate the cdf
    data = slices(:,ii:ii+2*window_length);
    
    %------------------------------------------------
    % Temp section to store data size
    
    data_size(ii,1) = ii;
    data_size(ii,2) = (getfield(whos('data'),'bytes'))/1024^2;    
    
    data = data(:);
    if min(data) ~= max(data)
        data = data(data ~= 0);                                 % remove hard zeros
        data = data(~isnan(data));                              % remove NaN's
        numdata = length(data);                                 % total number of real datapoints (without zeros and NaN's)
        widthbin = 3.5*std(data)/(numdata^(1/3));               % gives optimum bin width for gaussian distribution (fairly robust for other distributions)
        bins = ((min(data):widthbin:max(data)))';               % defines bin centres based on bin width and data extremes
        data_tmp = interp1(bins,bins,data,'nearest','extrap');  % interpolates data onto closest bin center
        check_bins = ~ismember(bins,data_tmp);                  % finds empty bins
        data_tmp = sort([data_tmp;bins(check_bins)]);           % temporarily pads data with empty bins to ensure no bins are empty, and sorts low to high
        [~,datapdf_tmp,~] = unique(data_tmp);                   % returns the highest index of the unique values in data
        datapdf_tmp = diff([0;datapdf_tmp]);                    % difference between indices is the number of values in a bin
        datapdf_tmp(check_bins) = 0;                            % reset the padded bins to zero
        datapdf = datapdf_tmp;                                  % put into a cell array
        datacdf = cumsum(datapdf)-datapdf/2;                    % cumulative sum across bins and re-centers to bin center
        datacdf = datacdf/max(datacdf);                         % normalises max of CDF to 1
        
        % Calculate probabilities from cdf
        data = slices(:,ii+window_length);                      % get central slice in window
        data = interp1(bins,bins,data,'nearest','extrap');      % interpolates data onto closest bin center
        [~,idxdataprob] = ismember(data,bins);                  % finds the bin index that the data is in
        
        start_ind = 1;
        end_ind = 0;
        for i_row = 1:1:size(trace_ilxl_bytes,1)            
            end_ind = end_ind + n_traces_row(i_row);
            dataprob(skeyn*(i_row-1)+live_traces(i_row,1):skeyn*(i_row-1)+live_traces(i_row,2),ii) ...
                = 1-datacdf(idxdataprob(start_ind:end_ind));      % lookup the probability scores from the cdf
            start_ind = end_ind + 1;
        end
    end
end

% Save the result to a binary file
fid = fopen(strcat(job_meta.output_dir,'sas_result_slice_',i_slab,'_slice.bin'),'w');
fwrite(fid,dataprob,'float32');
fclose(fid);

fid = fopen(strcat(job_meta.output_dir,'sas_result_slice_',i_slab,'_trace.bin'),'w');
fwrite(fid,dataprob','float32');
fclose(fid);

%% Gather data back to trace rather than slice order to write out to SEGY.

% Transpose input data
% dataprob_trans = dataprob';
% 
% results_out{1,1} = 'ilxl numbers';
% results_out{1,2}{1,1} = ilxl_read;
% results_out{2,1} = 'sas_result_trace';
% results_out{2,2} = dataprob_trans;
% 
% node_segy_write_traces(results_out,str2num(i_slab),job_meta.output_dir);

% write_blocks = 10000;
% % Split data into n_blocks to be written out.
% for kk = 1:1:write_blocks;
%     
%     % Calculate the trace to start reading from for this block
%     start_trace_idx = 1+((kk-1)*floor(n_traces/write_blocks);
%     
%     % Calculate the last trace to read for this block (adjusts for extra traces on last block)
%     if kk == n_blocks
%         end_trace_idx = n_traces;
%     else
%         end_trace_idx = kk*floor(n_traces/n_blocks);
%     end
%     
%     % Calculate the number of traces to read in this block
%     %n_traces_to_read = end_trace_idx-start_trace_idx+1;
%     
%     dataprob_trc = dataprob_trans(:,start_trace_idx:end_trace_idx);
%     fid = fopen(strcat(job_meta.output_dir,'sas_result_trace_',i_slab,'.bin'),'w');
%     fwrite(fid,dataprob_trc,'float32');  
% end
% fclose(fid);

end