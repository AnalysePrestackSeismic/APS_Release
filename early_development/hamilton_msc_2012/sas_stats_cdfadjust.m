function[average,variance] = sas_stats_cdfadjust(data_file_in,time_file_in,bytespersample,nt,ntwin,olap)
%
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
%
% PART 1.......................
%


%% READS IN AND CALCULATES SIZE OF DATA

% Open the input binary files.
fid1 = fopen(data_file_in);
fid2 = fopen(time_file_in);

% Work out the size of the binary files in bytes.
fseek(fid1, 0, 'eof'); % Moves to end of file
filesize1 = ftell(fid1); % filesize1 = currentposition(end of file)
frewind(fid1); %returns to the beginning of the file

%Repeats for input file 2
fseek(fid2, 0, 'eof');
filesize2 = ftell(fid2);
frewind(fid2);

%% Calculate the number of traces in the file
ntrace1 = filesize1/(nt*bytespersample);
% total number of samples in file/(number of time sample per trace * bytes per sample)
ntrace2 = filesize2/(nt*bytespersample);





%% READS IN THE DATA IN SLICE FORMAT FROM TIME VOLUME

% Read the first slice from the time volume. Function reads data in a
% binary format into column vector, the second input gives the size of the
% vector and the 3rd the precision, interprets based on the data format

tdud = fread(fid2,ntrace2,'float32');
t1 = fread(fid2,ntrace2,'float32');

% Read the second slice from the time volume.NOT SURE IF THIS CONTINUES
% FROM THE LAST OR IS THE SAME DATA
t2 = fread(fid2,ntrace2,'float32');



%% Initiates some simple variables.

skip_index = 1;
count = 1;


ntrace2

ntrace1


%% Splits data into slices which are 4ms apart, thus generating a 4ms sampling rate. Slices closer are note used. 4ms is the sampling rate in time and thus samples closer than this must have been interpolated and do not correspond
%% to real points at real time samples

while isempty(t2) ~= 1;  % while t2 i.e time slice2 is not empty isempty returns a 0 thus until it is empty the following sequence is initiated.   
    
    
  %  Corrects zeros and undefined values to be NaN.
        
    for k = 1:ntrace2                      % Cycles through every trace
     
        if t1(k) == 0 || t1(k) > 1e29           % tests the first trace in t1 to see if it 0 or really large
            t1(k) = NaN;                        % if so sets to NAN
        end
        
        if t2(k) == 0 || t2(k) > 1e29           % Repeats for t2
            t2(k) = NaN;
        end
        
        if t2(k) < t1(k)                        % Now compares t2 and t1 to determine is t2 is less than t1. if so then both t2 and t1 become NAN
            t2(k) = NaN;
            t1(k) = NaN;
        end
        
     end
    
    
%% Check to see if t2 and t1 are 4ms apart - required to have real sampling rate.     
% Check to see if t1 and t2 are at least 4ms apart.
% if the nan mean of t2 and t1 is greater than or equal to 4ms(nanmean
% gives the mean once NAN numbers have been removed.) 
%% THIS IF STATEMENT IS FAILING FOR SOME REASON
    
     if notnumbermean(t2) - notnumbermean(t1) >= 0.004
        
        % If t1 and t2 are at least 4ms apart then write to index the
        % number of samples between t1 and t2. first iteration = 0
        
        index(skip_index) = count-skip_index;
                
        % Set t1 to be equal to t2 ready for next loop.
        t1 = t2;
        
        % Read the next slice from the time volume ready for next loop.
        t2 = fread(fid2,ntrace2,'float32');
                
        % Set count equal to skip_index. This tracks the sample number of t2.
        count = skip_index;
                
        % Increment skip_index by one. This tracks the sample number of t1.
        skip_index=skip_index+1;
    
    else
        
        % If t1 and t2 are at not 4ms apart then read the next slice from the time volume ready for next loop.
        t2 = fread(fid2,ntrace2,'float32');
        
    end
    % Increment count by one. This tracks the sample number of t2.
    count = count+1;
end


%% Now the program has generated a vector with a range of indices which correspond to the locations of "true" timeslices within the volume. 

% Calculate the number of bytes needed to be skipped so that only slices which are at least 4ms apart are used in the average and variance calculations.
% index contains the number of time slices to jump this is translated into
%samples.



%Number of bytes to skip between reading in slices. 

skip = index*ntrace1*bytespersample;



%% Sample required for shift = size of time window(samples) - (size of window*overlapfraction) - in this case floor is used to round the window down. 


    
% Calculate the number of time samples to shift the calcualtion window by on each iteration such that olap is honoured.
win_shift = ntwin-(floor(ntwin*olap));

% Set the total number of loops needed. number of loops is number of skips
% - 
% is number of time slices. 
nloop = floor((length(skip)-ntwin)/win_shift);






%% PART 2..................

% Initialise some variables.
average = zeros(nloop+1,1);
variance = zeros(nloop+1,1);
stddev = zeros(nloop+1,1);

% Update the user on progress.
fprintf('\nCalculating statistics...\n');


%% Read the first window of data from data_file_in.
for i=1:ntwin
    data(((i-1)*ntrace1)+1:i*ntrace1,1) = fread(fid1,ntrace1,sprintf('%d*float32',ntrace1),skip(i));
end

length(data)
for i = 1:ntwin*ntrace1
    if data(i) > 1e29
        data(i) = NaN;
    end
end



% Calculate the statistics for the first window of data.
average(1) = notnumbermean(data);
variance(1) = notnumbervar(data); 
stddev(1) = sqrt(variance(1));



%% Loop to calculate statistics from each window of data.


for k = 1:nloop   
        clearvars data_tmp
        % If on the last loop then only read the remaining data from data_file_in.
    if k == nloop
        for i = 1:length(skip)-(ntwin+(win_shift*(k-1)))
            data_tmp(((i-1)*ntrace1)+1:i*ntrace1,1) = fread(fid1,ntrace1,sprintf('%d*float32',ntrace1),skip(ntwin+(win_shift*(k-1))+i));
        end
    % If not on the last loop read the extra data for the next time window.
    else
        for i = 1:win_shift
            data_tmp(((i-1)*ntrace1)+1:i*ntrace1,1) = fread(fid1,ntrace1,sprintf('%d*float32',ntrace1),skip(ntwin+(win_shift*(k-1))+i));
        end
    end
    % Correct undefined values to be NaN.    
    for i = 1:length(data_tmp)
        if data_tmp(i) > 1e29
            data_tmp(i) = NaN;
        end
    end
    
    % Combine previous data with that just read to produce the new time window of data.
    data = data(win_shift*ntrace1+1:end);
   
    % adds the newly read data onto the old data ( doing it this way simply
    % gives a more efficient way of reading in the data. 
    data = [data;data_tmp];

    % Calculate the statistics for the window of data.
    average(k+1) = notnumbermean(data);
    variance(k+1) = notnumbervar(data); 
    stddev(k+1) = sqrt(variance(k+1));
end



%% PART 3
% The purpose of this section is to give every point in the whole volume.
% Not only the time sampled slices an average and an std deviation.
% This section of the data gives average values to the whole data set. 
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
averagei = (interp1(x,average,xi,'pchip','extrap'));
stddevi = (interp1(x,stddev,xi,'pchip','extrap'));





%% PART 4
% Update the user on progress.
fprintf('\nCalculating cdf...\n');

% Initialise some variables.
data = zeros(ntrace1,1);


% Go back to the begining of data_file_in.
frewind(fid1);

% Create and open a file to output the anomaly indicator data.

fid3 = fopen(sprintf('cdfadjust1_%s',data_file_in),'a');

%ntwin,ntwin-win_shift,data_file_in),'a');




% Loop to calculate the anomaly indicator values and write them to file.



for i = 1:nt
    % Read a slice from data_file_in.
    data(:,1) = fread(fid1,ntrace1,'float32');
    
    [row1, col1] = find(data<averagei(i));
    [row2, col2] = find(data>=averagei(i));
   

%      if i == 3
%      row1
%      row2
%      length(row1)
%      length(row2)
%      end
%      mother = [row1; row2];
    
    
    prob1 = -1*(2*cdf1('Normal', abs(data(row1)-averagei(i)) + averagei(i),averagei(i),stddevi(i))-1);
    prob2 = (2*cdf1('Normal', abs(data(row2)-averagei(i)) + averagei(i),averagei(i),stddevi(i))-1);
   
    
    prob(row1,1) = prob1;
    prob(row2,1) = prob2;
    
    if i == 1
       
        length(prob)
        
    end
              
   %Assign every value on that slice a 'probability of being anomalous' based on the previous calculated statistics and assuming a Gaussian model.
   %prob = (2*cdf1('Normal',abs(data-averagei(i))+averagei(i),averagei(i),stddevi(i)))-1;
    
   %Write the anomally indicator value to file.
        
   fwrite(fid3,prob,'float32');
   
   %fwrite(fid3,prob2,'float32');
      
    %Update the user on progress.
        
    if single((nt/50)*floor(i/(nt/50))) == i;
            fprintf('%d / %d = %.0f%%\n',i,nt,(i*100)/(nt));
    end
end

% Close all files.
fclose all;

end



