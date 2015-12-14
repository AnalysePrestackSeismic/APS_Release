function wavelet_plot_norm_estimations_for_a_block(job_meta_path)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
plot_on = 0;
startvol = 1;
volinc = 1;
%endvol = job_meta.nvols;
endvol = 5;
totalvol = length(startvol:volinc:endvol);

% Load job meta information 
job_meta = load(job_meta_path);

%% Load wavelet set
wavelets = load(strcat(job_meta.wav_directory,'all_wavelets_time.mat'));
ns_wavelet = size(wavelets.all_wavelets_time{1},1)-1;
hns_wavelet = floor(ns_wavelet/2);
ns = job_meta.n_samples{1};
%i_block = str2double(i_block);
vol_count = 1;

for i_vol = startvol:volinc:endvol
    wavelet_z_grid = wavelets.all_wavelets_freq{i_vol}(1,:);
    wavelet = wavelets.all_wavelets_freq{i_vol}(2:end,:);
    if plot_on == 1
        figure(2)
        subplot(1,totalvol,vol_count); imagesc(wavelet);
        figure(3)
        subplot(1,totalvol,vol_count+1); plot(wavelet(:,10));
    end
    start_interp = min(wavelet_z_grid);%-hns_wavelet;
    end_interp = max(wavelet_z_grid);%+hns_wavelet;
    wavelet_interp{vol_count} = interpolate_wavelets(wavelet_z_grid,wavelet,start_interp,end_interp);
    wavelet_interp{vol_count} = circshift(ifft(wavelet_interp{vol_count}','symmetric'),floor(job_meta.ns_win/2));
    %wavelet_interp{i_vol} = wavelet_interp{i_vol};
    
    if start_interp > 1;
        % Pad with zeros or first wavelet
        pad_size = start_interp-1;
        wavelet_interp{vol_count} = [repmat(wavelet_interp{vol_count}(:,1),1,pad_size),...
            wavelet_interp{vol_count}];
    end
    if end_interp < ns;
        % Pad with zeros or last wavelet
        pad_size = ns-end_interp;
        wavelet_interp{vol_count} = [wavelet_interp{vol_count},...
            repmat(wavelet_interp{vol_count}(:,end),1,pad_size)];
    end
    if floor(ns_wavelet/2) == ns_wavelet/2
        wavelet_interp{vol_count}(end+1,:) = wavelet_interp{vol_count}(end,:);
    end
    if plot_on == 1
        figure(4)
        subplot(1,totalvol,vol_count); imagesc(wavelet_interp{vol_count});
    end
    vol_count = vol_count + 1;
end
interp_wavelet_z_grid = 1:ns;
ns_wavelet = size(wavelet_interp{1},1);
wavelet_norm = wavelet_rms_norm(totalvol,wavelet_interp,interp_wavelet_z_grid,ceil(totalvol*0.6667));
if plot_on == 1
    figure(5)
    imagesc(cell2mat(wavelet_norm))

end
%subplot(3,job_meta.nvols,i_vol+job_meta.nvols); plot(freq_axis,w_freq(3:2+ceil(job_meta.ns_win/2),:,i_vol)); set(gca,'XTick',1:ceil(10*freq_axis(2)):ceil(freq_axis(end)))

figure(10)

xmin = 1;
xmax = 129;
wavelets_vec = cell2mat(wavelet_norm);
wavelets_vec = wavelets_vec(:);
ymin = min(wavelets_vec);
ymax = max(wavelets_vec);


for ii = 1:1:ns
    
   subplot(1,totalvol,1); plot(wavelet_norm{1}(:,ii)); axis([xmin xmax ymin ymax])
   subplot(1,totalvol,2); plot(wavelet_norm{2}(:,ii)); axis([xmin xmax ymin ymax])
   subplot(1,totalvol,3); plot(wavelet_norm{3}(:,ii)); axis([xmin xmax ymin ymax])
   subplot(1,totalvol,4); plot(wavelet_norm{4}(:,ii)); axis([xmin xmax ymin ymax])
   subplot(1,totalvol,5); plot(wavelet_norm{5}(:,ii)); axis([xmin xmax ymin ymax])
   title(num2str(ii))
   pause(0.1)
end


% %for i_vol = 1:1:job_meta.nvols
%     fid_wav = fopen(strcat(job_meta.wav_directory,'fft_wavelets_block_',i_block,'.bin'));
%     w_tmp = fread(fid_wav,'float32');
%     fclose(fid_wav);
%     %n_win = w_tmp(1);
%     n_vol = w_tmp(1);
%     n_win = w_tmp(2);
%     %w_freq{i_vol} = reshape(w,[],n_win);
%     %w_freq{i_vol} = reshape(w_tmp(3:end),[],n_vol);
%     w_freq = reshape(w_tmp(3:end),[],n_win,n_vol);
%     % w_freq{i_vol} = w(3:end,:);
%     % Convert from frequency domain to time domain
%     w_time = circshift(ifft(w_freq(3:end,:,:),'symmetric'),floor(job_meta.ns_win/2));
% %end








% %imagesc(w_time{1});
% freq_axis = (1e6/job_meta.s_rate)/2*linspace(0,1,job_meta.ns_win/2);
% figure
% for i_vol = 1:1:job_meta.nvols
% subplot(3,job_meta.nvols,i_vol); imagesc(w_time(:,:,i_vol));
% subplot(3,job_meta.nvols,i_vol+job_meta.nvols); plot(freq_axis,w_freq(3:2+ceil(job_meta.ns_win/2),:,i_vol)); set(gca,'XTick',1:ceil(10*freq_axis(2)):ceil(freq_axis(end)))
% subplot(3,job_meta.nvols,(i_vol+(job_meta.nvols*2))); imagesc(w_freq(3:2+ceil(job_meta.ns_win/2),:,i_vol)');
% end
% %colorbar
% %colormap('gray')

end

function wavelet_interp = interpolate_wavelets(wavelet_z_grid,wavelet,start_interp,end_interp)
%     start_interp = min(wavelet_z_grid)-mode(diff(wavelet_z_grid'));
%     end_interp = max(wavelet_z_grid)+mode(diff(wavelet_z_grid'));
    wavelet_interp = interp1(wavelet_z_grid,wavelet',start_interp:1:end_interp,'linear','extrap');
end

function wavelet_norm = wavelet_rms_norm(totalvol,wavelet_interp,wavelet_z_grid,vol_index)
    % Normalise the wavelets to have constant energy w.r.t. angle. The energy
    % is set to that of the nearest angle wavelets. Wavelet energy still varies
    % w.r.t. time.
    A = cell2mat(wavelet_interp);
    %B = sqrt(sum(A.^2));
    B = sqrt(max(A.^2));
    C = reshape(B',length(wavelet_z_grid),[]);
    D = C(:,vol_index);
    for ii=1:totalvol
        E = A(:,1+(ii-1)*length(wavelet_z_grid):ii*length(wavelet_z_grid));
        F = bsxfun(@rdivide,bsxfun(@times,E,D'),sqrt(sum(E.^2)));
        wavelet_norm{ii} = F;
    end
end





