function [] = seismic_anomaly_spotter_lite(seismic_mat_path,i_block,n_blocks,window_length,binary_slice_path,output_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    i_block = str2double(i_block);
    n_blocks = str2double(n_blocks);
    window_length = str2double(window_length);

    % Get some information about the seimsic data
    tmp_seismic = fread(fopen(seismic_mat_path,'r'),'double');
    seismic.n_samples = tmp_seismic(1,1);
    seismic.trace_ilxl_bytes = reshape(tmp_seismic(2:end),3,[])';
    seismic.n_traces = size(seismic.trace_ilxl_bytes,1);
    
    % Work out the min/max ils/xls for the survey bounding box
    ilmin = min(seismic.trace_ilxl_bytes(:,1));
    ilmax = max(seismic.trace_ilxl_bytes(:,1));
    ilinc = mode(diff(unique(seismic.trace_ilxl_bytes(:,1))));
    iln = 1+((ilmax-ilmin)/ilinc);
    
    xlmin = min(seismic.trace_ilxl_bytes(:,2));
    xlmax = max(seismic.trace_ilxl_bytes(:,2));
    xlinc = mode(diff(unique(seismic.trace_ilxl_bytes(:,2))));
    xln = 1+((xlmax-xlmin)/xlinc);
    
    % Make ils/xls of the survey bounding box
    bounding_box = [reshape(repmat((ilmin:ilinc:ilmax),xln,1),[],1),repmat((xlmin:xlinc:xlmax)',iln,1),(1:1:iln*xln)'];
    
    % Get the linear indices of the lives traces in the survey bounding box
    live_traces = bounding_box(ismember(bounding_box(:,1:2),seismic.trace_ilxl_bytes(:,1:2),'rows'),3);
     
    % Read the slices for this slab
    [slices] = node_binary_read_slices_lite(seismic,i_block,n_blocks,binary_slice_path,window_length);
    
    % Work out the number of windows to loop through
    n_windows = size(slices,2)-2*window_length;
    
    % Pre-allocate memory for the output zero-padded data probability slab
    dataprob = zeros(iln*xln,n_windows);

    for ii = 1:n_windows
        % Calculate the cdf
        data = slices(:,ii:ii+2*window_length);
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
            dataprob(live_traces,ii) = 1-datacdf(idxdataprob);      % lookup the probability scores from the cdf
        end
    end
        
    % Save the result to a binary file
    fid = fopen(strcat(output_dir,sprintf('sas_result_slices_block_%d.bin',i_block)),'w');
    fwrite(fid,dataprob,'float32');
    fclose(fid);
    
    % Write trace order blocks (number of blocks to divide output into)
    
end