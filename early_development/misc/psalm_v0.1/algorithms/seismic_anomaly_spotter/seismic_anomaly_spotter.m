function [] = seismic_anomaly_spotter(block_mat_all,block_mat,process_files_mat,output_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    % process_positions = node_make_processing_positions(block_mat);
    % load(process_files_mat);
    process_files = load(process_files_mat);
    [traces, process_positions] = segy_read_multiple_files(block_mat_all,process_files_mat, block_mat);
    [~, process_positions_ilxl_index] = ismember(process_positions.ilxl_grid,process_positions.ilxl_pos,'rows');
    % [~, process_positions_z_index] = ismember(process_positions.z_grid,process_positions.z_pos,'rows');
    for ii = 1:length(process_positions.z_grid)
        if ii == 1
            [~, process_positions_z_index_tmp] = ismember(process_positions.z_grid(ii,1),flipud(process_positions.z_pos),'rows');
        else
            [~, process_positions_z_index_tmp] = ismember(process_positions.z_grid(ii,1),flipud(process_positions.z_pos(process_positions_z_index(ii-1,1)+process_positions.z_step+1:end,1)),'rows');
        end
        process_positions_z_index(ii,1) = length(process_positions.z_pos)-process_positions_z_index_tmp+1; 
    end

    for ii = 1:process_files.nfiles
        traces_tmp = traces{ii};
        for kk = 1:length(process_positions_z_index)
            if process_positions_z_index(kk)-process_positions.z_aperture <= 0
                first_z = 1;
            else
                first_z = process_positions_z_index(kk)-process_positions.z_aperture;
            end
            if process_positions_z_index(kk)+process_positions.z_aperture > length(process_positions.z_pos)
                last_z = length(process_positions.z_pos);
            else
                last_z = process_positions_z_index(kk)+process_positions.z_aperture;
            end
            subvol_traces_tmp = traces_tmp(first_z:last_z,:)'; %window of traces to be used below

            % Calculate the cdf
            data = subvol_traces_tmp(:);
            if min(data) == max(data)
                dataprob(kk,:) = 0;   
            else
                data = data(data ~= 0);                                         % remove hard zeros
                data = data(~isnan(data));                                      % remove NaN's
                numdata = length(data);                                         % total number of real datapoints (without zeros and NaN's)
                widthbin(kk,1) = 3.5*std(data)/(numdata^(1/3));                 % gives optimum bin width for gaussian distribution (fairly robust for other distributions)
                bins(:,kk) = {(single(min(data):widthbin(kk,1):max(data)))'};   % defines bin centres based on bin width and data extremes
                data_tmp = interp1(bins{kk},bins{kk},data,'nearest','extrap');  % interpolates data onto closest bin center
                check_bins = ~ismember(bins{kk},data_tmp);                      % finds empty bins
                data_tmp = sort([data_tmp;bins{kk}(check_bins)]);               % temporarily pads data with empty bins to ensure no bins are empty, and sorts low to high
                [~,datapdf_tmp,~] = unique(data_tmp);                           % returns the highest index of the unique values in data
                datapdf_tmp = diff([0;datapdf_tmp]);                            % difference between indices is the number of values in a bin
                datapdf_tmp(check_bins) = 0;                                    % reset the padded bins to zero
                datapdf(:,kk) = {datapdf_tmp};                                  % put into a cell array
                datacdf(:,kk) = {cumsum(datapdf{kk})-datapdf{kk}/2};            % cumulative sum across bins and re-centers to bin center 
                datacdf(:,kk) = {single(datacdf{kk}/max(datacdf{kk}))};         % normalises max of CDF to 1
%                 average(kk,1) = mean(data(~isnan(data)));                      % calculates average of the slab

                data = subvol_traces_tmp(:,1+round((last_z-first_z)/2)); %finds central slice in window
%                 data = data-average(kk);
                data = interp1(bins{kk},bins{kk},data,'nearest','extrap');
                [~,idxdataprob] = ismember(data,bins{kk});
                idxdataprob = single(idxdataprob);
                dataprob(kk,:) = datacdf{kk}(idxdataprob)';
    %             dataprob(kk,dataprob(kk,:)<0.5) = 1-dataprob(kk,dataprob(kk,:)<0.5);
    %             dataprob(kk,:) = (2*dataprob(kk,:))-1;
                dataprob(kk,:) = 1-dataprob(kk,:);
            end
        end
    end
    results_out{1,1} = 'dataprob';
    results_out{1,2} = double(dataprob);
    
    save(strcat(output_dir,sprintf('%s_results_block_%d',process_files.func_name,process_positions.block_id)),'results_out','-v7.3');
end