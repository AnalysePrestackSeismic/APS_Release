%function [] = dynava(slices, slice_file_dir, ref_angle, wave_length, I_file_out, G_file_out)
%

warning off

loc = pwd;
cd(slice_file_dir);
    
% Figure out the number of files in the current directory
[~,nfiles] = (system('ls -B | wc -l')); 
nfiles = str2double(nfiles);

% Read all filenames and convert from ascii to double
[~,fnames] = system('ls -B1');
numeric = double(fnames);

% Preallocate memory for some variables
count = 2;
fname_index = zeros(1,nfiles+1);

% Loop to separate out each file name from one long character string
for ij= 1:length(fnames)
    if numeric(1,ij) == 10
        fname_index(1,count) = ij;
        count = count+1;
    end
end

% Loop to read each file as ascii into a cell array
for ik=1:nfiles
    files_in{ik} = fnames(1,(fname_index(1,ik)+1):(fname_index(1,ik+1)-1));
end

fprintf('\nThe following files have been found:\n')
for il = 1:nfiles
    fprintf('File %d: %s\n',il,files_in{il});
    fid_in(il,1) = fopen(files_in{il});
end

cd(loc);

% Open the input binary files.
fid_I = fopen(I_file_out,'a');
fid_G = fopen(G_file_out,'a');

win = ((3*wave_length)/4)+1;
stepout = (win-1)/2;

for i = 1+stepout:slices{1}.n_samples-stepout
    
    % Calculate the wavelets
    for k = 1:nfiles
        fseek(fid_in(k),slices{k}.slice_pointers(i-stepout,2),'bof');
        data{k} = fread(fid_in(k),[slices{k}.n_traces,win],'*float32');
        data{k}(data{k} > 1e29) = NaN;
        wavelet(:,k) = circshift(ifft(mean(abs(fft(data{k}')),2)),floor(win/2));
        if slices{k}.angle == ref_angle
            ref = k;
        end
    end
    
    % Calculate and apply matching filters
    for k = 1:nfiles        
        matchfilt(:,k) = match_call(wavelet(:,ref),wavelet(:,k),-6);
        data{k} = conv2(data{k}',matchfilt(:,k),'same');
        data{k} = data{k}';
        angles(k,1) = slices{k}.angle;
    end
    
    % Fit the data using 2-term AVA
    ava = nan(2,slices{1}.n_traces);
    
    if i == 1+stepout
        for l = 1:1+stepout
            for k = 1:nfiles
                data_ava(:,k) = data{k}(:,l);
            end
            ava = [ones(size(angles)) sin(angles*pi/180)]\data_ava';
            fwrite(fid_I,ava(1,:),'float32');
            fwrite(fid_G,ava(2,:),'float32');
        end
    elseif i == slices{1}.n_samples-stepout
        for l = slices{1}.n_samples-stepout:slices{1}.n_samples
            for k = 1:nfiles
                data_ava(:,k) = data{k}(:,l-slices{1}.n_samples+win);
            end
            ava = [ones(size(angles)) sin(angles*pi/180)]\data_ava';
            fwrite(fid_I,ava(1,:),'float32');
            fwrite(fid_G,ava(2,:),'float32');
        end
    else
        for k = 1:nfiles
            data_ava(:,k) = data{k}(:,stepout+1);
        end
        ava = [ones(size(angles)) sin(angles*pi/180)]\data_ava';
        fwrite(fid_I,ava(1,:),'float32');
        fwrite(fid_G,ava(2,:),'float32');
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

%end

