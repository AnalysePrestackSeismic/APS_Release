function[] = sas_stats(data_file_in,time_file_in,bytespersample,nt,ntwin,olap)
% Calculates slice-consistent windowed statistics for a slice-optimised seismic volume (as output by the function optslice).
%
% USAGE:   [average,variance] = sas_stats(data_file_in,time_file_in,bytespersample,nt,ntwin,olap)
%
% INPUTS:   data_file_in = slice-optimised binary volume (as output by the function optslice) of attribute data.
%           time_file_in = slice-optimised binary volume (as output by the function optslice) of time data in seconds.
%           bytespersample = number of bytes per sample (4 if input data is float32 format).
%           nt = number of time samples in input data.
%           ntwin = calculation window length in number of time samples.
%           olap = fractional overlap of sliding time windows (0 = none, 0.1 = 10%, etc).
%
% OUTPUTS:  average = vector of average values within each window of data from data_file_in.
%           variance = vector of variance values within each window of data from data_file_in.
%           cdf_win_(ntwin)_olap_(ntwin-winshift)_(data_file_in) = binary file of slice optimised anomaly indicator volume. This can be made to inline optimised using the function optil.

% Open the input binary files.
fid1 = fopen(data_file_in);
fid2 = fopen(time_file_in);

% Work out the size of the binary files in bytes.
fseek(fid1, 0, 'eof');
filesize1 = ftell(fid1);
frewind(fid1);

fseek(fid2, 0, 'eof');
filesize2 = ftell(fid2);
frewind(fid2);

% Calculate the number of traces in the files.
ntrace1 = filesize1/(nt*bytespersample);
ntrace2 = filesize2/(nt*bytespersample);

% Read the first slice from the time volume.
t1 = fread(fid2,ntrace2,'*float32');

% Read the second slice from the time volume.
t2 = fread(fid2,ntrace2,'*float32');

% Initialise some variables.
skip_index = 1;
count = 1;

% Begin loop to evaluate which slices are at least 4ms apart. Slices that are closer than this are not used to calculate statistics because their data must have been interpolated. Run loop until the lower slice (t2) is empty.
while isempty(t2) ~= 1;
    % Corrects zeros and undefined values to be NaN.
    t1(t1 == 0 | t1 > 1e29 | t2 < t1) = NaN;
    t2(t2 == 0 | t2 > 1e29 | t2 < t1) = NaN;
%     for k = 1:ntrace2
%         if t1(k) == 0 || t1(k) > 1e29
%             t1(k) = NaN;
%         end
%         if t2(k) == 0 || t2(k) > 1e29
%             t2(k) = NaN;
%         end
%         if t2(k) < t1(k)
%             t2(k) = NaN;
%             t1(k) = NaN;
%         end
%     end
    % Check to see if t1 and t2 are at least 4ms apart.
    if mean(t2(~isnan(t2))) - mean(t1(~isnan(t1))) >= 0.004
        % If t1 and t2 are at least 4ms apart then write to index the number of samples between t1 and t2.
        index(skip_index) = count-skip_index;
        % Set t1 to be equal to t2 ready for next loop.
        t1 = t2;
        % Read the next slice from the time volume ready for next loop.
        t2 = fread(fid2,ntrace2,'*float32');
        % Set count equal to skip_index. This tracks the sample number of t2.
        count = skip_index;
        % Increment skip_index by one. This tracks the sample number of t1.
        skip_index=skip_index+1;
    else
        % If t1 and t2 are at not 4ms apart then read the next slice from the time volume ready for next loop.
        t2 = fread(fid2,ntrace2,'*float32');
    end
    % Increment count by one. This tracks the sample number of t2.
    count = count+1;
end

% Calculate the number of bytes needed to be skipped so that only slices which are at least 4ms apart are used in the avaergae and variance calculations.
skip = index*ntrace1*bytespersample;

% Calculate the number of time samples to shift the calcualtion window by on each iteration such that olap is honoured.
win_shift = ntwin-(floor(ntwin*olap));
% Set the total number of loops needed.
nloop = floor((length(skip)-ntwin)/win_shift);

% Initialise some variables.
average = zeros(nloop+1,1);
%variance = zeros(nloop+1,1);
%stddev = zeros(nloop+1,1);

% Update the user on progress.
fprintf('\nCalculating statistics...\n');

% Read the first window of data from data_file_in.
for i=1:ntwin
    data(((i-1)*ntrace1)+1:i*ntrace1,1) = fread(fid1,ntrace1,sprintf('%d*float32=>float32',ntrace1),skip(i));
end

% Correct undefined values to be NaN.
data(data > 1e29) = NaN;
% for i = 1:ntwin*ntrace1
%     if data(i) > 1e29
%         data(i) = NaN;
%     end
% end

% Build pdf and cdf from data (assumes data symmetry about zero)
idxnotnan = ~isnan(data);
numdata = sum(idxnotnan);
%numbin = ceil(sqrt(numdata));
%numbin = 1+2*ceil(numbin/2);
%widthbin(1,1) = (max(data)-min(data))/numbin;
widthbin(1,1) = 3.5*std(data(idxnotnan))/(numdata^(1/3));
%bins{1,1} = widthbin(1,1)*(0.5-numbin/2):widthbin(1,1):widthbin(1,1)*(-0.5+numbin/2);
bins{1,1} = single(min(data):widthbin(1,1):max(data));
datapdf{1,1} = hist(data,bins{1,1});
datacdf{1,1} = cumsum(datapdf{1,1})-datapdf{1,1}/2;  
%lengthdatacdf(1,1) = length(datacdf{1,1});
%maxdatacdf(1,1) = max(datacdf{1,1});
datacdf{1,1} = single(datacdf{1,1}/max(datacdf{1,1}));

% Calculate the statistics for the first window of data.
average(1) = mean(data(~isnan(data)));
%variance(1) = var(data(~isnan(data))); 
%stddev(1) = sqrt(variance(1));

% Loop to calculate statistics from each window of data.
for k = 1:nloop    
    clearvars data_tmp
    % If on the last loop then only read the remaining data from data_file_in.
    if k == nloop
        for i = 1:length(skip)-(ntwin+(win_shift*(k-1)))
            data_tmp(((i-1)*ntrace1)+1:i*ntrace1,1) = fread(fid1,ntrace1,sprintf('%d*float32=>float32',ntrace1),skip(ntwin+(win_shift*(k-1))+i));
        end
    % If not on the last loop read the extra data for the next time window.
    else
        for i = 1:win_shift
            data_tmp(((i-1)*ntrace1)+1:i*ntrace1,1) = fread(fid1,ntrace1,sprintf('%d*float32=>float32',ntrace1),skip(ntwin+(win_shift*(k-1))+i));
        end
    end
    % Correct undefined values to be NaN.
    data_tmp(data_tmp > 1e29) = NaN;
%     for i = 1:length(data_tmp)
%         if data_tmp(i) > 1e29
%             data_tmp(i) = NaN;
%         end
%     end
    
    % Combine previous data with that just read to produce the new time window of data.
    data = data(win_shift*ntrace1+1:end);
    data = [data;data_tmp];
    
    % Build pdf and cdf from data (assumes data symmetry about zero)
    idxnotnan = ~isnan(data);
    numdata = sum(idxnotnan);
    %numbin = ceil(sqrt(numdata));
    %numbin = 1+2*ceil(numbin/2);
    %widthbin(k+1,1) = (max(data)-min(data))/numbin;
    widthbin(k+1,1) = 3.5*std(data(idxnotnan))/(numdata^(1/3));
    %bins{k+1,1} = widthbin(k+1,1)*(0.5-numbin/2):widthbin(k+1,1):widthbin(k+1,1)*(-0.5+numbin/2);
    bins{k+1,1} = single(min(data):widthbin(k+1,1):max(data));
    datapdf{k+1,1} = hist(data,bins{k+1,1});
    datacdf{k+1,1} = cumsum(datapdf{k+1,1})-datapdf{k+1,1}/2;
    %lengthdatacdf(k+1,1) = length(datacdf{k+1,1});
    %maxdatacdf(k+1,1) = max(datacdf{k+1,1});
    datacdf{k+1,1} = single(datacdf{k+1,1}/max(datacdf{k+1,1}));

    % Calculate the statistics for the window of data.
    average(k+1) = mean(data(~isnan(data)));
    %variance(k+1) = var(data(~isnan(data))); 
    %stddev(k+1) = sqrt(variance(k+1));
end

% Pad data cdfs and bins with nans so they are all the same length.
% regdatacdf = nan(nloop+1,max(lengthdatacdf));
% regbins = nan(nloop+1,max(lengthdatacdf));
% idxcdf = (max(lengthdatacdf) - lengthdatacdf)/2;
% for i = 1:nloop+1
%     regdatacdf(i,idxcdf(i)+1:max(lengthdatacdf)-idxcdf(i)) = datacdf{i,1};
%     regdatacdf(i,:) = regdatacdf(i,:)/max(regdatacdf(i,:));
%     regbins(i,idxcdf(i)+1:max(lengthdatacdf)-idxcdf(i)) = bins{i,1};
% end

% regdatacdf = regdatacdf';
% regbins = regbins';

% Update the user on progress.
fprintf('\nInterpolating...\n');

% Work out the sample numbers for the data used in the statistics calculations.
cumu(1) = 1;
for i=2:length(index)+1
    cumu(i)=cumu(i-1)+(index(i-1)+1);
end

% Interpolate the statistics onto the full data time range.
x = downsample(cumu(ntwin/2:end-ntwin/2),win_shift);
xi = 1:1:nt;
averagei = single(interp1(x,average,xi,'pchip','extrap'));
%stddevi = (interp1(x,stddev,xi,'pchip','extrap'));
%datacdfi = (interp1(x,regdatacdf,xi,'nearest','extrap'));
%binsi = (interp1(x,regbins,xi,'nearest','extrap'));
datacdfloc = interp1(x,x,xi,'nearest','extrap');
datacdfidx = diff(datacdfloc)./diff(datacdfloc);
datacdfidx(isnan(datacdfidx)) = 0;
datacdfidx = [0,datacdfidx];
datacdfidx = 1+cumsum(datacdfidx);
datacdfidx = logical(datacdfidx');

% Update the user on progress.
fprintf('\nCalculating cdf...\n');

clearvars -except fid1 ntwin win_shift data_file_in ntrace1 datacdfidx datacdf bins averagei fid4 nt

% Initialise some variables.
%data = zeros(ntrace1,1);

% Go back to the begining of data_file_in.
frewind(fid1);

% Create and open a file to output the anomaly indicator data.
%fid3 = fopen(sprintf('modelprob_win%d_olap%d_%s',ntwin,ntwin-win_shift,data_file_in),'a');
fid4 = fopen(sprintf('dataprob_win%d_olap%d_%s',ntwin,ntwin-win_shift,data_file_in),'a');

% Loop to calculate the anomaly indicator values and write them to file.
for i = 1:nt
    % Read a slice from data_file_in.
    data(:,1) = fread(fid1,ntrace1,'*float32');
    % Assign every value on that slice a 'probability of being anomalous' based on the previous calculated statistics and assuming a Gaussian model.
%    modelprob = (2*cdf('Normal',abs(data-averagei(i))+averagei(i),averagei(i),stddevi(i)))-1;
    %[~,idxdataprob] = min(abs(bsxfun(@minus,bins{datacdfidx(i)}',data')));
    data = abs(data-averagei(i))+averagei(i);
    [~,idxdataprob] = min(abs(bsxfun(@minus,bins{datacdfidx(i)}',data')));
    %dataprob = datacdf{datacdfidx(i)}(idxdataprob);
    dataprob = (2*(datacdf{datacdfidx(i)}(idxdataprob)))-1;
    % Write the anomally indicator value to file.
%    fwrite(fid3,modelprob,'float32');
    fwrite(fid4,dataprob,'float32');
    % Update the user on progress.
    if single((nt/50)*floor(i/(nt/50))) == i;
            fprintf('%d / %d = %.0f%%\n',i,nt,(i*100)/(nt));
    end
end

% Close all files.
fclose all;

end