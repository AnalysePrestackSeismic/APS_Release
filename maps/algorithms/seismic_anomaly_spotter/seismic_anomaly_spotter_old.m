function [] = seismic_anomaly_spotter_old(job_meta_path,i_slab,vol_index,n_slab,window_length)
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

pkey_inc_mode = mode(job_meta.pkey_inc);
skey_inc_mode = mode(job_meta.skey_inc);

pkeyn = 1+((job_meta.pkey_max(str2double(vol_index))-job_meta.pkey_min(str2double(vol_index)))...
    /job_meta.pkey_inc(str2double(vol_index)));
skeyn = 1+((job_meta.skey_max(str2double(vol_index))-job_meta.skey_min(str2double(vol_index)))...
    /job_meta.skey_inc(str2double(vol_index)));

start_slice_idx = (1+((str2double(i_slab)-1)*floor(job_meta.n_samples{str2double(vol_index)}...
    /str2double(n_slab))))-str2double(window_length);

% Calculate the last trace to read for this slab (adjusts for extra slices in last slab)
if str2double(i_slab) == str2double(n_slab)
    end_slice_idx = str2double(window_length)+job_meta.n_samples{str2double(vol_index)};
else
    end_slice_idx = str2double(window_length)+(str2double(i_slab)*floor(job_meta.n_samples{str2double(vol_index)}/str2double(n_slab)));
end
pad_start = 0;
pad_end = 0;
% vol_traces needs to be padded size
% Correct for window overshoots
if start_slice_idx <= 0
    pad_start = 1+abs(start_slice_idx);
    %pad_start = zeros(seismic.n_traces/n_blocks,1+abs(start_slice_idx));
    start_slice_idx = 1;
end
if end_slice_idx >= job_meta.n_samples{str2double(vol_index)}
    % if end_slice_idx >= ns
    pad_end = end_slice_idx-job_meta.n_samples{str2double(vol_index)};
    end_slice_idx = job_meta.n_samples{str2double(vol_index)};
    %     end_slice_idx = ns;
end

% Calculate the number of slices to read in this slab
n_slices_to_read = end_slice_idx-start_slice_idx+1;
% Pre-allocate memory for the output zero-padded data probability slab
% generate volume for storing slices
slices = zeros(pkeyn*skeyn,pad_start+n_slices_to_read+pad_end);
% win_idx = (0:1:str2double(window_length)-1)';
% Read traces for each block
% Extract time samples and populate slices slab
for i_block = 1:1:job_meta.n_blocks
    [~, traces, ilxl_read, ~] = ...
        node_segy_read(job_meta_path,vol_index,num2str(i_block));
    if size(traces) > 1
        % Loop round to pick the water bottom        
        %wb_idx = water_bottom_picker(traces,0);
        %wb_idx(wb_idx < 0) = 1;
        wb_idx = ones(1,size(traces,2));
        win_sub = bsxfun(@plus,wb_idx,(start_slice_idx-1:end_slice_idx-1)');
        win_sub_log = win_sub<job_meta.n_samples{str2double(vol_index)};
        win_sub = win_sub(1:min(sum(win_sub_log)),:);
        win_ind = bsxfun(@plus,win_sub,(0:job_meta.n_samples{str2double(vol_index)}:...
           job_meta.n_samples{str2double(vol_index)}*(size(traces,2)-1)));
        
        %win_ind = win_ind(:,sum(win_sub > job_meta.n_samples{str2double(vol_index)})==0);
        if isempty(win_ind)
            fprintf('Window indices exceed trace length\n')
            return
        end
%        figure
%        imagesc(traces(:,1:1000))
%        hold all
%        plot(win_sub(1,1:1000))
%        plot(win_sub(end,1:1000))
       
        %traces = trace_flatten(traces,wb_idx{i_block});
        % Position these as linear indices
        %fprintf('%d\n', i_block)
        n_iline = (ilxl_read(:,1)-job_meta.pkey_min(str2double(vol_index)))/pkey_inc_mode+1;
        n_xline = (ilxl_read(:,2)-job_meta.skey_min(str2double(vol_index)))/skey_inc_mode+1;
        lin_ind = ((n_iline-1).*skeyn)+n_xline;

        % Trim traces
        traces = traces(win_ind);
        n_slices_to_read = size(traces,1);
        slices(double(lin_ind),pad_start+1:pad_start+n_slices_to_read) = traces';
    end
end
clearvars traces ilxl_read n_iline n_xline pad_start pad_end lin_ind

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
        data = interp1(bins,bins,data,'nearest','extrap');  % interpolates data onto closest bin center
        check_bins = ~ismember(bins,data);                  % finds empty bins
        data = sort([data;bins(check_bins)]);           % temporarily pads data with empty bins to ensure no bins are empty, and sorts low to high
        [~,datapdf,~] = unique(data);                   % returns the highest index of the unique values in data
        datapdf = diff([0;datapdf]);                    % difference between indices is the number of values in a bin
        datapdf(check_bins) = 0;                            % reset the padded bins to zero
        %datapdf = datapdf_tmp;                                  % put into a cell array
        datacdf = cumsum(datapdf)-datapdf/2;                    % cumulative sum across bins and re-centers to bin center
        datacdf = datacdf/max(datacdf);                         % normalises max of CDF to 1
        
        % Calculate probabilities from cdf
                            % get central slice in window
        data = interp1(bins,bins,slices(:,ii+window_length),'nearest','extrap');      % interpolates data onto closest bin center
        [~,idxdataprob] = ismember(data,bins);                  % finds the bin index that the data is in
        
        dataprob(:,ii) = 1-datacdf(idxdataprob);  
    end
end

% Save the result to a binary file
% fid = fopen(strcat(job_meta.output_dir,'sas_result_slice_',i_slab,'_slice.bin'),'w');
% fwrite(fid,dataprob,'float32');
% fclose(fid);

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