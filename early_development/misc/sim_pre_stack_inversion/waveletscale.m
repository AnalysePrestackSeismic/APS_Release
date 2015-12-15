function [scaled_wavelet_data] = waveletscale(input_angle_stacks_dir, log_data, wavelet_data, seisrange, c1, c2, c3)

start_point = pwd;
cd(input_angle_stacks_dir)

nwells = length(log_data);
[~,nwavelets] = size(wavelet_data);
nwavelets = nwavelets-1;

% Figure out the number of files in the current directory
[~,nfiles] = (system('ls -B | wc -l')); 
nfiles=str2double(nfiles);

if nfiles ~= nwavelets
    error('There must be the same number of wavelets as there are angle stack volumes. Currently there are %d wavelets but %d angle stack volumes!',nwavelets,nfiles);
end

% Read all filenames and convert from ascii to double
[~,fnames] = system('ls -B1');
numeric = double(fnames);

% Preallocate memory for some variables
count = 2;
fname_index = zeros(1,nfiles+1);
trace_data{nwells,1} = [];
syn_data{nwells,1} = [];
ref_data{nwells,1} = [];
scalar = zeros(nwells,nfiles);
scaled_wavelet_data = zeros(size(wavelet_data));

% Loop to separate out each file name from one long character string
for i = 1:length(fnames)
    if numeric(1,i) == 10
        fname_index(1,count) = i;
        count = count+1;
    end
end

% Loop to read traces from each angle stack at each well location as binary
% float32 into a cell array and to create synthetic seismograms for each
% angle
for i = 1:nwells
    
    % Calculate number of bytes to skip having corrected for log times not
    % falling on seismic samples
    traceskip = ((log_data{i,1}(1,1)-seisrange(1,1))/seisrange(1,3))*(((seisrange(2,2)-seisrange(2,1))/seisrange(2,3))+1);
    timeskip = ((seisrange(3,2)-seisrange(3,1))/seisrange(3,3))+1;
    extraskip = (((log_data{i,1}(1,2)-seisrange(2,1))/seisrange(2,3))*(timeskip))+((log_data{i,1}(1,3)-seisrange(3,1))/seisrange(3,3));
    tread = round(log_data{i,1}(1,4)/seisrange(3,3))-round(log_data{i,1}(1,3)/seisrange(3,3))+2;
    byteskip = 4*((traceskip*timeskip)+extraskip);
    
    ref_imp = zeros((length(log_data{i,1}(:,1)))-1,3);
    for l = 2:4
        ref_imp(:,l-1) = [diff(log_data{i,1}(2:end,l));0]./conv2(log_data{i,1}(2:end,l),[1,1],'same');   
    end
    ref_imp = ref_imp(1:end-1,:);
    
    for k = 1:nfiles
        
        % Work out angle stack filenames
        file_in = fnames(1,(fname_index(1,k)+1):(fname_index(1,k+1)-1));
        angle = str2double(file_in(1,1:2));
        if isnan(angle)
            error('Angle stack filenames must begin with the angle of the stack e.g ''05_stack'' or ''25_stack''!')
        end
        
        % Read trace at well location over log time range from angle stacks
        fid = fopen(file_in);
        fseek(fid,byteskip,'bof');
        trace_data{i,1}(:,k) = [angle;fread(fid,tread,'float32')];
        frewind(fid);
        fclose(fid);
        
        % Calculate angle dependent reflectivity and accompanying synthetic
        % seismograms
        ref_data{i,1}(:,k) = [angle;(c1(k,1)*ref_imp(:,1))+(c2(k,1)*ref_imp(:,2))+(c3(k,1)*ref_imp(:,3))];
        syn_data{i,1}(:,k) = [angle;conv2(ref_data{i,1}(2:end,k),wavelet_data(2:end,k+1),'same')];
    end
    trace_data{i,1} = sortrows(trace_data{i,1}',1)';
end
% this loop was not in correct place!

% Calculate scalars
for k = 1:nfiles
    scalar(i,k) = sqrt(sum(trace_data{i,1}(2:end,k).^2))/(sqrt(sum(ref_data{i,1}(2:end,k).^2))*sqrt(sum(wavelet_data(2:end,k+1).^2)));
end
avg_scalar = mean(scalar,1);
scaled_wavelet_data(:,1) = wavelet_data(:,1);
scaled_wavelet_data(1,:) = wavelet_data(1,:);
for k = 1:nfiles
    scaled_wavelet_data(2:end,k+1) = wavelet_data(2:end,k+1)*avg_scalar(1,k);
end

% Go back to the directory you were in to start with
cd(start_point);
end

