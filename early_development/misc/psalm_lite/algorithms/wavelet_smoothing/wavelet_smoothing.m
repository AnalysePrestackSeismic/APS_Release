function [] = wavelet_smoothing(wavelet_mat_path)
%WAVELET_SMOOTHING average the wavelets over some space

wavlim = 26;

% find all the Wavelets in a direcotry ###################################

strcomm = ['ls -B1  ',wavelet_mat_path,'average_wavelets_time*.mat | wc -l'];
strcommb = ['ls -B1  ',wavelet_mat_path,'average_wavelets_time*.mat'];
[~,nfiles] = (system(strcomm));
nfiles=str2double(nfiles);

[~,fnames] = system(strcommb);
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
    wavelet_mat_pathtime = fnames(1,(fname_index(1,ik)+1):(fname_index(1,ik+1)-1));
    % files_in.path{cumm_files+ik-1} = input_dir{ii};
    load(wavelet_mat_pathtime,'avg_w_time');
    all_wavelets_time{1,ik} = avg_w_time;
    clear avg_w_time;
end
plotwav(nfiles,1);
% maybe l;imit the output range
for aj = 1:nfiles 
    all_wavelets_time{1,aj} = all_wavelets_time{1,aj}(1:end,1:wavlim);
end
save(strcat(char(wavelet_mat_path),'all_wavelets_time'),'all_wavelets_time','-v7.3');

% read in the freq domain wavelets ################################

strcomm = ['ls -B1  ',wavelet_mat_path,'average_wavelets_freq*.mat | wc -l'];
strcommb = ['ls -B1  ',wavelet_mat_path,'average_wavelets_freq*.mat'];
[~,nfiles] = (system(strcomm));
nfiles=str2double(nfiles);

[~,fnames] = system(strcommb);
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
    wavelet_mat_pathtime = fnames(1,(fname_index(1,ik)+1):(fname_index(1,ik+1)-1));
    % files_in.path{cumm_files+ik-1} = input_dir{ii};
    load(wavelet_mat_pathtime,'avg_w');
    all_wavelets_freq{1,ik} = avg_w;
    clear avg_w;
end


% maybe l;imit the output range
for aj = 1:nfiles 
    all_wavelets_freq{1,aj} = all_wavelets_freq{1,aj}(1:end,1:wavlim);
end
save(strcat(char(wavelet_mat_path),'all_wavelets_freq'),'all_wavelets_freq','-v7.3');


plotwav(nfiles,2);

% Plot Wavelets ###################################################
    function [] = plotwav(ns,figno)
    % parameters are number of plots and figure number        

        figure(figno)
        
        for a = 1:ns  
            subplot(8,1,a);
            imagesc(bsxfun(@rdivide,all_wavelets_time{1,a}(2:end,:),max(all_wavelets_time{1,a}(2:end,:))));
        end
        
        print 
        
        
    end

end

