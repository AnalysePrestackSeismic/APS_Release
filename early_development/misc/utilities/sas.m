function [] = sas(slices, slice_file_in, stepout, dataprob_file_out)
%Calculates cdf probabilities from amplitudes in slice_file_in. A cdf is
%created for the data within a window of length 2*stepout+1 about each
%slice and used to assign probabilities to the amplitudes on the current 
%slice.

% Open the input binary files.
fid_in = fopen(slice_file_in);
fid_out = fopen(dataprob_file_out,'a');
win = (2*stepout)+1;

for i = 1+stepout:slices{1}.n_samples-stepout
    % Calculate the cdf
    fseek(fid_in,slices{1}.slice_pointers(i-stepout,2),'bof');
    data = fread(fid_in,[slices{1}.n_traces*win,1],'*float32');
    data(abs(data) > 1e10) = NaN;
    idxnotnan = ~isnan(data);
    numdata = sum(idxnotnan);
    widthbin(i-stepout:i+stepout,1) = 3.5*std(data(idxnotnan))/(numdata^(1/3));
    bins(i-stepout:i+stepout,1) = {single(min(data):widthbin(1,1):max(data))};
    data_tmp = sort(interp1(bins{i,1},bins{i,1},data,'nearest','extrap'));
    check_bins = ~ismember(bins{i,1},data_tmp);
    data_tmp = sort([data_tmp;bins{i,1}(check_bins)']);
    [~,datapdf_tmp,~] = unique(data_tmp);
    datapdf_tmp = diff([0;datapdf_tmp]);
    datapdf_tmp(check_bins) = 0;
    datapdf(i-stepout:i+stepout,1) = {datapdf_tmp};
    datacdf(i-stepout:i+stepout,1) = {cumsum(datapdf{i,1})-datapdf{i,1}/2};  
    datacdf(i-stepout:i+stepout,1) = {single(datacdf{i,1}/max(datacdf{i,1}))};
%     average(i-stepout:i+stepout,1) = mean(data(~isnan(data)));

    % Classify the data
    if i == 1+stepout
        data = data(1:(stepout+1)*slices{1}.n_traces,1);
%         data = abs(data-average(i))+average(i);
        data = interp1(bins{i,1},bins{i,1},data,'nearest','extrap');
        data(isnan(data)) = max(bins{i,1});
        [~,idxdataprob] = ismember(data,bins{i,1});
        idxdataprob = single(idxdataprob);
%         dataprob = ((2*(datacdf{i,1}(idxdataprob)))-1)';
        dataprob = 1-datacdf{i,1}(idxdataprob)';
        fwrite(fid_out,dataprob,'float32');
    elseif i == slices{1}.n_samples-stepout
        data = data(1+(stepout*slices{1}.n_traces):end,1);
%         data = abs(data-average(i))+average(i);
        data = interp1(bins{i,1},bins{i,1},data,'nearest','extrap');
        data(isnan(data)) = max(bins{i,1});
        [~,idxdataprob] = ismember(data,bins{i,1});
        idxdataprob = single(idxdataprob);
%         dataprob = ((2*(datacdf{i,1}(idxdataprob)))-1)';
        dataprob = 1-datacdf{i,1}(idxdataprob)';
        fwrite(fid_out,dataprob,'float32');
    else
        data = data(1+(stepout*slices{1}.n_traces):(stepout+1)*slices{1}.n_traces,1);
%         data = abs(data-average(i))+average(i);
        data = interp1(bins{i,1},bins{i,1},data,'nearest','extrap');
        data(isnan(data)) = max(bins{i,1});
        [~,idxdataprob] = ismember(data,bins{i,1});
        idxdataprob = single(idxdataprob);
%         dataprob = ((2*(datacdf{i,1}(idxdataprob)))-1)';
        dataprob = 1-datacdf{i,1}(idxdataprob)';
        fwrite(fid_out,dataprob,'float32');
    end

    % Progress update
    if floor(i/10) == i/10;
            fprintf('Processing %d of %d (%d%% complete)\n',i,slices{1}.n_samples,round((i/slices{1}.n_samples)*100));
    end
end

fid_geom = fopen('geom.bin','a');
fwrite(fid_geom,slices{1}.geometry(1,1),'float32');
fwrite(fid_geom,slices{1}.geometry(2,1),'float32');
fwrite(fid_geom,slices{1}.geometry(3:end,1),'int');

fclose all;

end

