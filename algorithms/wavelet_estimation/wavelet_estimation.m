function [] = wavelet_estimation(job_meta_path,i_block,decimate,first_live_block)
%% ------------------ Disclaimer  ------------------
% 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) makes no representation or warranty, express or implied, in 
% respect to the quality, accuracy or usefulness of this repository. The code
% is this repository is supplied with the explicit understanding and 
% agreement of recipient that any action taken or expenditure made by 
% recipient based on its examination, evaluation, interpretation or use is 
% at its own risk and responsibility.
% 
% No representation or warranty, express or implied, is or will be made in 
% relation to the accuracy or completeness of the information in this 
% repository and no responsibility or liability is or will be accepted by 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) in relation to it.
%% ------------------ License  ------------------ 
% GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
%% github
% https://github.com/AnalysePrestackSeismic/
%% ------------------ FUNCTION DEFINITION ---------------------------------
% wavelet_estimation: Function to estimate wavelets for DIGI
% Pick water bottom on 2/3 of max angle stack to get reliable pick
% Use this pick to flatten all other angle stacks
% Estimate wavelets from flattened angle stacks
%   Arguments:
%       job_meta_path:      Path to job meta file as a string
%       i_block:            The block on which the program will be run
%       decimate:           decimation factor
%       first_live_block:   The block number of the first live block
%
%   Outputs:
%	Saves to Disk: Wavelet file for ith block in binary format
%   Writes to Disk:
%       nothing

plot_on = 0;
decimate = str2double(decimate);
% Load job meta information
lsout = ls(job_meta_path);
job_meta = load(job_meta_path);

% Wavelet estimation parameters
ns_win = 128;       % Number of Samples in the wavelet
%ns_overlap = 32;
ns_overlap = 96;    % Number of Samples overlap between moving wavelet window
padding = 50;


% Make directory to save results
if exist(strcat(job_meta.output_dir,'wavelets/'),'dir') == 0
    mkdir(strcat(job_meta.output_dir,'wavelets/'))
end
wav_directory = strcat(job_meta.output_dir,'wavelets/');

% this bit needs to be changes, might as well stack all the
% offsets/angles/angle stacks and then pick wb on those, this would make
% both post and pre stack reads the same and put the data into a 3d array
% Read traces for 2/3 max angle stack
%vol_index_wb = ceil(job_meta.nvols*0.6667);
if job_meta.is_gather == 0
    pick_wb_ind = ceil(job_meta.nvols*0.6667);
    vol_index_wb = 1;
    [~, traces{vol_index_wb}, ilxl_read{vol_index_wb}] = node_segy_read(job_meta_path,num2str(pick_wb_ind),i_block);
    % check to make sure it read something if not exit
    if size(traces{vol_index_wb},1) == 1 && size(traces{vol_index_wb},2) == 1
        return
    end
else
    % find the water bottom on a few gathers and take the middle one t use
    % as an offset plane to pick the water bottom on
    vol_index_wb = 1;
    % read all the data for this block
    % node_segy_read(job_meta_path,vol_index,i_block)
    [~, gathers, ilxl_read, offset_read] = node_segy_read(job_meta_path,'1',i_block);
    % find the total number of offsets
    offset = unique(offset_read);
    if isempty(gathers) == 1 &&  isempty(ilxl_read) == 1 && isempty(offset_read) == 1
        return
    end
    %     juststack = 1;
    %     if juststack == 1
     traces{vol_index_wb} = zeros(size(gathers,1),size(gathers,2)/size(offset,2)); 
    for stki  = 1:size(offset,2)
        traces{vol_index_wb} = traces{vol_index_wb} + gathers(:,offset_read == offset(stki));
    end
    % should reshape and sum instead of this terrible loop
    %     else
    % this part of the code could be used to make a sub stack around
    % the wb pick time that is the interpolate to find the peak for an
    % accurate wb picker
    %read the middle 5 gathers in the input data to find the middle tkey(angle/offset) value of the water bottom
    %traces{vol_index_wb} = gathers(:,(floor(size(gathers,2)/2)-(length(offset)*2)):(floor(size(gathers,2)/2)+(length(offset)*2)));
    %traces{vol_index_wb} = gathers(:,(floor(size(gathers,2)/2)-(length(offset)*2.5)+2):(floor(size(gathers,2)/2)+(length(offset)*2.5)));
    %tmpoffread = offset_read((floor(size(gathers,2)/2)-(length(offset)*2.5)+2):(floor(size(gathers,2)/2)+(length(offset)*2.5)));
    traces{vol_index_wb} = gathers(:,floor((size(gathers,2)/2)-(length(offset)*2.5)+2):floor((size(gathers,2)/2)+(length(offset)*2.5)));
    tmpoffread = offset_read(floor((size(gathers,2)/2)-(length(offset)*2.5)+2):floor((size(gathers,2)/2)+(length(offset)*2.5)));
    %pick the water bottom
    [wb_idxcj] = water_bottom_picker(traces{vol_index_wb}(:,:),0);
    %filter the water bottom pick to make a difference in WB time for
    %traces that are not picking the wb
    filtw = [1 2 3 2 1]/9;
    wb_idxcjfilt = conv(wb_idxcj,filtw,'same');
    % now make a blank array the size of the wb index array
    wb_idxcj2 = zeros(1,size(wb_idxcj,2));
    %now find the difference between the values of the filtered and
    %unfiltered water bottom indexes, if there is a difference then it is
    %likely to not be the water bottom as it should be mostly flat on the
    %gathers
    wb_idx_diff_ck = abs((wb_idxcj./wb_idxcjfilt)-1);
    wb_idxcj2(wb_idx_diff_ck < 0.009) =  wb_idxcj(wb_idx_diff_ck < 0.009);
    
    %now work out the index locations of the water bottom
    wb_idx_index = 1:1:size(wb_idxcj,2);
    %apply a logical index to the index array to give the index of where the wb is picked and less then 10 elsewhere
    wb_idx_index(ismember(wb_idxcj2,floor((min(wb_idxcj2((wb_idxcj2 > 10)))*0.9)):1:ceil((min(wb_idxcj2((wb_idxcj2 > 10)))*1.1))));
    % select the angles with the wb on
    tmpwbangs = (tmpoffread(wb_idx_index(ismember(wb_idxcj2,floor((min(wb_idxcj2((wb_idxcj2 > 10)))*0.9)):1:ceil((min(wb_idxcj2((wb_idxcj2 > 10)))*1.1))))));
    tmpangpickstd = ceil(std(double(tmpwbangs)));
    tmpangpick = floor(mean(tmpwbangs));
    %tmpangpick = floor(mean(tmpoffread(wb_idx_index(ismember(wb_idxcj2,floor((min(wb_idxcj2((wb_idxcj2 > 10)))*0.9)):1:ceil((min(wb_idxcj2((wb_idxcj2 > 10)))*1.1)))))));
    pick_wb_ind = find(offset == tmpangpick);
    %tmp_first = max([min(offset) (tmpangpick - tmpangpickstd)]);
    %tmp_last = min([max(offset) (tmpangpick + tmpangpickstd)]);
    %         tmp_first_wbidx = find(offset == max([min(offset) (tmpangpick - tmpangpickstd)]));
    %         tmp_last_wbidx = find(offset == min([max(offset) (tmpangpick + tmpangpickstd)]));
    %         %ismember(offset_read,(tmpangpick - tmpangpickstd):(tmpangpick + tmpangpickstd)
    %
    %         traces{vol_index_wb} = zeros(size(gathers,1),size(gathers,2)/size(offset,2));
    %         for stki  = tmp_first_wbidx:tmp_last_wbidx
    %             traces{vol_index_wb} = traces{vol_index_wb} + gathers(:,offset_read == offset(stki));
    %         end
    %traces{vol_index_wb} = traces{vol_index_wb}./(tmp_last_wbidx - tmp_first_wbidx + 1);
    
    %     end
    %     wb_idx_index(wb_idxcj2 == min(wb_idxcj2((wb_idxcj2 > 10))));
    %     % calculate the gather index in each gather by removing the integer
    %     % number of gathers from the index number and putting back to all being
    %     % the same angle index in each gather ie 1-46,1-46,1-46 etc... and
    %     % findingf the average value of all the indexes
    %     pick_wb_ind = floor(mean(((wb_idx_index(wb_idxcj2 == min(wb_idxcj2((wb_idxcj2 > 10)))))/length(offset) - floor((wb_idx_index(wb_idxcj2 == min(wb_idxcj2((wb_idxcj2 > 10)))))/length(offset)))*length(offset)));
    %     %read the offset plane from the input data gathers.
    traces{vol_index_wb} = gathers(:,offset_read == offset(pick_wb_ind));
end

%

% [~, traces, ~,~] = ...
%     node_segy_read(job_meta_path,num2str(vol_index_wb),i_block);
if plot_on == 1 && size(traces{1},2) > 1
    figure(2)
    imagesc(traces{1}(:,:));
end
% Pick water bottom

if isfield(job_meta, 'wb_path')
    %wb_idx = dlmread(job_meta.wb_path,'delimiter','\t');
    wb_idx_in = dlmread(job_meta.wb_path);
    % col 1 inline
    % col 2 xline
    % col 3 twt
    [~,locations] = ismember(ilxl_read{1}(1:end,:),wb_idx_in(:,1:2),'rows');
    %wb_idx = zeros(size(traces{vol_index_wb},2),1);
    zero_loc = locations ~= 0;
    %wb_idx(zero_loc) = wb_idx_in(locations(zero_loc),3);
    xi = (1:size(traces{vol_index_wb},2))';
    x = xi(zero_loc)';
    wb_idx = interp1(x,wb_idx_in(locations(zero_loc),3),xi);
    clear wb_idx_in
    %min_il = min(ilxl_read{vol_index_wb}(:,1));
    %max_il = max(ilxl_read{vol_index_wb}(:,1));
    %min_xl = min(ilxl_read{vol_index_wb}(:,2));
    %max_xl = max(ilxl_read{vol_index_wb}(:,2));
    %wb_idx2 = wb_idx(wb_idx(:,1) >= min_il & wb_idx(:,1) <= max_il & wb_idx(:,2) >= min_xl & wb_idx(:,2) <= max_xl,:);
    %    wb_idx = wb_idx(wb_idx(:,1) >= min_il & wb_idx(:,1) <= (max_il+1) & wb_idx(:,2) >= min_xl & wb_idx(:,2) <= (max_xl+1),:);
    %   wb_idx(:,3) = wb_idx(:,3)./(job_meta.s_rate/1000);
    wb_idx = (wb_idx./(job_meta.s_rate/1000))';
    %padding = 10;
    wb_idx = round(wb_idx-padding);
    win_sub = bsxfun(@plus,wb_idx,(0:job_meta.n_samples{vol_index_wb}-max(wb_idx))');
    
    win_ind = bsxfun(@plus,win_sub,(0:job_meta.n_samples{vol_index_wb}:...
        job_meta.n_samples{vol_index_wb}*(size(traces{vol_index_wb},2)-1)));
else
    [wb_idx] = water_bottom_picker(traces{vol_index_wb}(:,1:decimate:end),padding);
    wb_idx(wb_idx < 0) = 1;
    %wb_idx(isnan(wb_idx)) = 1;
    if job_meta.is_gather == 1
        %wbpickfile = fopen(strcat(wav_directory,'wbpick_',i_block,'.xyz'),'w');
        ilxltoprint = ilxl_read(1:length(offset):end,:);
        %dlmwrite(strcat(wav_directory,'wbpick_',i_block,'.xyz'),[ilxltoprint,(int32(wb_idx+10)'.*job_meta.s_rate)/1000],'delimiter', ' ', 'precision', '%-6d','newline', 'unix');
        % only write out the values that are not 1, so are picked
        dlmwrite(strcat(wav_directory,'wbpick_',i_block,'.xyz'),[ilxltoprint((wb_idx ~= 1)',:),(uint32(wb_idx(wb_idx ~= 1)+padding)'.*job_meta.s_rate)/1000],'delimiter', ' ', 'precision', '%-6d','newline', 'unix');
        %fprintf(wbpickfile,'%6d %6d %6d\n',ilxltoprint[,wb_idx);
        %fclose(wbpickfile);
        
    else
        ilxltoprint = ilxl_read{vol_index_wb};
        dlmwrite(strcat(wav_directory,'wbpick_',i_block,'.xyz'),[ilxltoprint((wb_idx ~= 1)',:),(uint32(wb_idx(wb_idx ~= 1)+padding)'.*job_meta.s_rate)/1000],'delimiter', ' ', 'precision', '%-6d','newline', 'unix');
    end
    win_sub = bsxfun(@plus,wb_idx,(0:job_meta.n_samples{vol_index_wb}-max(wb_idx))');
    win_ind = bsxfun(@plus,win_sub,(0:job_meta.n_samples{vol_index_wb}:...
        job_meta.n_samples{vol_index_wb}*(size(traces{vol_index_wb}(:,1:decimate:end),2)-1)));
end

% calculate water bottom index
wb_ind_avg = mean(wb_idx);


% Loop over all volumes and windows to estimate wavelets, so for angle
% stacks these are seperate volumes, for angle gathers they are different
% planes in the gathers volume

if job_meta.is_gather == 0
    i_vol_max = job_meta.nvols;
else
    i_vol_max = length(offset);
end

for i_vol = 1:1:i_vol_max
    if job_meta.is_gather == 0
        [~, traces, ~, ~] = node_segy_read(job_meta_path,num2str(i_vol),i_block);
    else
        traces = gathers(:,offset_read == offset(i_vol));
    end
    %[n_samples,n_traces] = size(traces);
    traces = traces(:,1:decimate:end);
    [~,n_traces] = size(traces);
    if n_traces > 2
        traces = traces(win_ind);
        %variance(i_vol) = var(traces(:));
        
        % updated taper to be the same as int_grad_inv
        %taper = linspace(0,1,25)';
        taper = single([ zeros(floor(padding-5),1)',linspace(0,1,20)])';
        %traces(1:25,:) = bsxfun(@times,traces(1:25,:),taper);
        traces(1:length(taper),:) = bsxfun(@times,traces(1:length(taper),:),taper);
        
        variance(i_vol) = var(traces(:));
        % cjedit_30122013
        [n_samples,n_traces] = size(traces); % n samples needs to be redefined as the data has been flatterned on the waterbottom index and the array truncated
        %
        % keep record of number of samples included in std
        % n_pop_std(i_vol) = n_samples*n_traces;
        start_index = 1:ns_win-ns_overlap-1:n_samples-ns_win;
        end_index = start_index+ns_win-1;
        %     if end_index(end) < n_samples
        %         end_index(end+1) = n_samples;
        %         start_index(end+1) = n_samples-ns_win+1;
        %     end
        n_win = length(start_index);
        w = zeros(2+ns_win,n_win);
        %
        % make taper to apply to signal before fft
        taperlen = 16;
        %taperst = linspace(0,1,taperlen)';
        taperst = (sin(linspace((-pi/2),(pi/2),taperlen)')+1)/2;
        taperend = 1 - taperst;
        taperapply = [taperst;ones((ns_win-(taperlen*2)),1);taperend];
        %
        for ii = 1:n_win
            % Estimate wavelets and store meta information
            %if ii == n_win
            %    w(1,ii) = n_samples;
            %else
            w(1,ii) = start_index(ii);
            %end
            %w(2,ii) = n_traces;
            %NFFT = 2^nextpow2(size(win_sub,1));
            % make a tmp array with the fft values in it
            tmpfft = zeros(2,2);
            tmpfft = abs(fft(bsxfun(@times,traces(start_index(ii):end_index(ii),:),taperapply)));
            %w(2,ii) = floor((sum(tmpfft(:) ~= 0,'double'))/ns_win);
            
            w(2,ii) = sum(sum(tmpfft ~= 0,1,'double') ~= 0,2);
            w(3:end,ii) = sum(tmpfft,2,'double');
            
            %w(3:end,ii) = sum(abs(fft(bsxfun(@times,traces(start_index(ii):end_index(ii),:),taperapply))),2);
            % index as linear index, summed not averaged!
            
            % Produce plots showing the windows
            if ii == 1 && plot_on == 1
                figure(1)
                imagesc(traces(:,1:1000))
                hold on
                plot(repmat(end_index(ii),1000),'--')
                hold off
            elseif plot_on == 1
                hold on
                plot(repmat(end_index(ii),1000),'--')
                hold off
                %figure(2)
                %subplot(n_win,1,ii);
                %imagesc(traces(start_index(ii):end_index(ii),:));
            end
            
        end
        %end
        % cjedit_30122013 - added if to not write out empty wavelets as far
        % too many files othwerwise
        %cjedit 02/01/2014 - commentout and made a single output file per
        %iblock below
        %    %if n_traces > 1
        %        % Save estimated wavelets for this volume
        %        fid_wav = fopen(strcat(wav_directory,job_meta.volumes{i_vol},'_fft_wavelets_block_',i_block,'.bin'),'w');
        %        fwrite(fid_wav,n_win,'float32');
        %        fwrite(fid_wav,w,'float32');
        %        fclose(fid_wav);
        %end
        % cjedit 02/01/2014 added a write to only one file per iblock and made
        % a slot for angle volume number in the output and total number of
        % angle volumes
        %%
        % Save estimated wavelets for this volume
        if i_vol == 1
            fid_wav = fopen(strcat(wav_directory,'fft_wavelets_block_',i_block,'.bin'),'w');
            %fwrite(fid_wav,job_meta.nvols,'float32');
            fwrite(fid_wav,i_vol_max,'float32');
            fwrite(fid_wav,n_win,'float32');
            fwrite(fid_wav,variance(1),'float32');  % this is used by wsmooth parameter in the inversion
            fwrite(fid_wav,pick_wb_ind,'float32');  % offset on which the water bottom was picked
            fwrite(fid_wav,wb_ind_avg,'float32');   % z of water bottom in block
        end
        %fwrite(fid_wav,i_vol,'float32');
        %fwrite(fid_wav,n_win,'float32');
        fwrite(fid_wav,w,'float32');
    end
    
end
if n_traces > 2
    fclose(fid_wav);
    % Open file again to write statistics
    fid_wav = fopen(strcat(wav_directory,'fft_wavelets_block_',i_block,'.bin'),'r');
    
    % Calculate
    stdev = round(sqrt(median(variance))); % might like to think about calculating this properly based on live samples
    w = fread(fid_wav,'float32');
    fclose(fid_wav);
    w(3) = stdev;
    fid_wav = fopen(strcat(wav_directory,'fft_wavelets_block_',i_block,'.bin'),'w');
    fwrite(fid_wav,w,'float32');
    fclose(fid_wav);
end

%If this is the first live block do the following write into job meta file
if str2double(i_block) == str2double(first_live_block)
    % Add processing information to job meta
    job_meta.wav_directory = wav_directory;
    job_meta.ns_win = ns_win;
    %job_meta.n_win = n_win; not the same for each wavelet file
    job_meta.ns_overlap = ns_overlap;
    save(job_meta_path,'-struct','job_meta','-v7.3');
end

end
