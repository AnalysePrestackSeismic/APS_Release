function [] = wavelet_estimation(block_mat_all,block_mat,process_files_mat,output_dir)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    maxNumCompThreads(1)
    process_files = load(process_files_mat);
    [traces, process_positions] = segy_read_multiple_files(block_mat_all,process_files_mat, block_mat); 
    
    for ii = 1:size(process_positions.ilxl_grid,1)
        if ii == 1
            [~, process_positions_ilxl_index_tmp] = ismember(process_positions.ilxl_grid(ii,1:2),flipud(process_positions.ilxl_pos(:,1:2)),'rows');
        else
            [~, process_positions_ilxl_index_tmp] = ismember(process_positions.ilxl_grid(ii,1:2),flipud(process_positions.ilxl_pos(process_positions_ilxl_index(ii-1,1)+process_positions.ilxl_step+1:end,1:2)),'rows');
        end
        process_positions_ilxl_index(ii,1) = length(process_positions.ilxl_pos)-process_positions_ilxl_index_tmp+1; 
    end
    for ii = 1:size(process_positions.z_grid,1)
        if ii == 1
            [~, process_positions_z_index_tmp] = ismember(process_positions.z_grid(ii,1),flipud(process_positions.z_pos),'rows');
        else
            [~, process_positions_z_index_tmp] = ismember(process_positions.z_grid(ii,1),flipud(process_positions.z_pos(process_positions_z_index(ii-1,1)+process_positions.z_step+1:end,1)),'rows');
        end
        process_positions_z_index(ii,1) = length(process_positions.z_pos)-process_positions_z_index_tmp+1; 
    end
    wavelet_length = 75;
    half_wavelet_length = floor(wavelet_length/2);
    time_sampling = 0.004;

    for ii = 1:process_files.nfiles
        traces_tmp = traces{ii};
        [nt cols] = size(traces_tmp);
        if process_positions.z_aperture == 0 % if no z aperture then use all z-samples
            if process_positions.ilxl_aperture == 0 % if no ilxl aperture then use all traces
                subvol_traces_tmp = traces_tmp;
                [subvol_traces_tmp, n_live_cols] = find_dead_traces(subvol_traces_tmp);
                ft_traces_tmp = (sum(abs(fft(subvol_traces_tmp,nt)),2))./n_live_cols;
                coords = [process_positions.ilxl_pos(process_positions_ilxl_index,1:2)';-1*(ones(1,length(process_positions_ilxl_index)))];
            else % if there is a ilxl aperture then use sub-selection of traces
                for kk = 1:length(process_positions_ilxl_index)
                    subvol_traces_tmp = traces_tmp(...
                        :,...
                        process_positions_ilxl_index(kk)-((2*(process_positions.ilxl_aperture^2))+2*process_positions.ilxl_aperture):...
                        process_positions_ilxl_index(kk)+((2*(process_positions.ilxl_aperture^2))+2*process_positions.ilxl_aperture));
                    [subvol_traces_tmp, n_live_cols] = find_dead_traces(subvol_traces_tmp);
                    ft_traces_tmp(:,kk) = (sum(abs(fft(subvol_traces_tmp,nt)),2))./n_live_cols;
                end
                coords = [process_positions.ilxl_pos(process_positions_ilxl_index,1:2)';-1*(ones(1,length(process_positions_ilxl_index)))];
            end
        else % if there is a z aperture then use sub-selection of z-samples
            nt = 1+2*process_positions.z_aperture;
            if process_positions.ilxl_aperture == 0  % if no ilxl aperture then use all traces
                % coords = repmat(process_positions.ilxl_pos(process_positions_ilxl_index,1:2)',1,length(process_positions_z_index));
                coords = [repmat(process_positions.ilxl_pos(round(length(process_positions.ilxl_pos)/2),1:2)',1,length(process_positions_z_index));process_positions_z_index'];
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
                    subvol_traces_tmp = traces_tmp(first_z:last_z,:);
                    [subvol_traces_tmp, n_live_cols] = find_dead_traces(subvol_traces_tmp);
                    ft_traces_tmp(:,kk) = (sum(abs(fft(subvol_traces_tmp,nt)),2))./n_live_cols;
                    % coords(3,1+(kk-1)*length(process_positions_ilxl_index):kk*length(process_positions_ilxl_index)) = process_positions_z_index(kk);
                end
            else % if there is a ilxl aperture then use sub-selection of traces
                coords = repmat(process_positions.ilxl_pos(process_positions_ilxl_index,1:2)',1,length(process_positions_z_index));
                for jj = 1:length(process_positions_z_index) 
                    for kk = 1:length(process_positions_ilxl_index)
                        if process_positions_z_index(jj)-process_positions.z_aperture <= 0
                            first_z = 1;
                        else
                            first_z = process_positions_z_index(jj)-process_positions.z_aperture;
                        end
                        if process_positions_z_index(jj)+process_positions.z_aperture > length(process_positions.z_pos)
                            last_z = length(process_positions.z_pos);
                        else
                            last_z = process_positions_z_index(jj)+process_positions.z_aperture;
                        end
                        subvol_traces_tmp = traces_tmp(first_z:last_z,...
                            process_positions_ilxl_index(kk)-((2*(process_positions.ilxl_aperture^2))+2*process_positions.ilxl_aperture):...
                            process_positions_ilxl_index(kk)+((2*(process_positions.ilxl_aperture^2))+2*process_positions.ilxl_aperture));
                        [subvol_traces_tmp, n_live_cols] = find_dead_traces(subvol_traces_tmp);
                        ft_traces_tmp(:,kk+(jj-1)*length(process_positions_ilxl_index)) =...
                            (sum(abs(fft(subvol_traces_tmp,nt)),2))./n_live_cols;
                    end
                    coords(3,1+(jj-1)*length(process_positions_ilxl_index):jj*length(process_positions_ilxl_index)) = process_positions_z_index(jj);
                end
            end
        end
        wavelet{1,ii} = make_wavelet(ft_traces_tmp,wavelet_length);
        wavelet{1,ii} = [coords;wavelet{1,ii}];
        clearvars ft_traces_tmp
    end
    time_axis = [(-half_wavelet_length*time_sampling:time_sampling:-time_sampling)';0;(time_sampling:time_sampling:half_wavelet_length*time_sampling)'];
    save(strcat(output_dir,'wavelet_block_',num2str(process_positions.block_id)),'wavelet','-v7.3');
end

%%
function [subvol_traces_tmp, n_live_cols] = find_dead_traces(subvol_traces_tmp)
    cols = size(subvol_traces_tmp,2);
    live_samples = ~isnan(subvol_traces_tmp);
    dead_traces_tmp = sum(live_samples);
    dead_traces = zeros(1,cols);
    dead_traces(dead_traces_tmp == 0) = 1;
    n_dead_traces = sum(dead_traces,2);
    subvol_traces_tmp(~live_samples) = 0;
    n_live_cols = (cols-n_dead_traces);
    n_live_cols(n_live_cols == 0) = nan;
end

%%
function [wavelet] = make_wavelet(ft_traces_tmp, wavelet_length)
    ft_taper_length = 10; % assymetric (tapers only the high frequency energy to zero at Nyquist frequency)
    wavelet_taper_length = 30; % symmetric about zero (outside +/- wavelet taper length about zero, amplitude is set to zero, while inside there is a cosine taper)
    half_wavelet_length = floor(wavelet_length/2);
    nt = size(ft_traces_tmp,1);
    hnt = floor(nt/2);
    ft_taper = [ones(hnt-ft_taper_length,1); ((1+cos((0:1/(ft_taper_length-1):1)*pi)')/2)];
    ft_taper = [ft_taper; zeros(nt-length(ft_taper),1)];
    ft_traces_tmp = bsxfun(@times,ft_traces_tmp,ft_taper); % ensures zero energy at Nyquist frequency
    wavelet = circshift(ifft(ft_traces_tmp,'symmetric'),hnt);
    [~, peak_index] = max(wavelet);
    peak_index = max(peak_index);
    wavelet = wavelet(peak_index-half_wavelet_length:peak_index+half_wavelet_length,:);
    wavelet_taper = [((1+cos((-1:1/(wavelet_taper_length-1):0)*pi)')/2); ones(13,1); ((1+cos((0:1/(wavelet_taper_length-1):1)*pi)')/2)];
    wavelet_taper = [zeros((size(wavelet,1)-length(wavelet_taper))/2,1); wavelet_taper; zeros((size(wavelet,1)-length(wavelet_taper))/2,1)];
    wavelet = bsxfun(@times,wavelet,wavelet_taper); % ensures zero amplitude at terminations
    wavelet(isnan(wavelet)) = 0;
end