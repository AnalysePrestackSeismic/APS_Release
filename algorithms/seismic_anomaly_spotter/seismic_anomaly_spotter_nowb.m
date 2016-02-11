function [] = seismic_anomaly_spotter(job_meta_path,vol_index,window_length,amp,hilb,flip_phase,start_slab,end_slab)
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
%%-------------------------------------------------------------------------
% Function Defination:
% SAS function creats an anolamlous amplitude probabality
% volume. It takes a volume and operates on it on a moving window basis 
% where the window length is defined by the user and the window is moved by
% unity after every computation step. In a window it finds the anomalous amplitudes...


% Input:
%   job_meta_path: Path of the job meta file
%   vol_index:  Volume index for which  volume to use
%   window_length: WIndow length in samples for the algorithm. 
%   amp: 
%   start_slab:Start sample number of the slab?
%   end_slab: End Sample number of the slab? (fail safe value)

% Output:
%   function has void output but it writes out
%   a volume of seismic anolmaly volume probablity 
%% Load Job Meta Information and Create the Framework
% Load job meta information 
job_meta = load(job_meta_path);                                 % Load job meta file

% Slab estimation parameters
ns_win = str2double(window_length);                             % convert window length to number
start_index = 1:1:job_meta.n_samples{str2double(vol_index)};    % Create array of Start sample index for all the windows
end_index = start_index+ns_win;                                 % Create array of End sample index for all the windows

pkey_inc_mode = mode(job_meta.pkey_inc);                        % Primary Key (inline) increment (mode )
skey_inc_mode = mode(job_meta.skey_inc);                        % Secondary Key (inline) Increment (mode)

pkeyn = 1+((job_meta.pkey_max(str2double(vol_index))-job_meta.pkey_min(str2double(vol_index)))...
    /job_meta.pkey_inc(str2double(vol_index)));                 % Calculate Number of inlines (primary key)
skeyn = 1+((job_meta.skey_max(str2double(vol_index))-job_meta.skey_min(str2double(vol_index)))...
    /job_meta.skey_inc(str2double(vol_index)));                 % Calculate Number of inlines (secondary key)

if str2double(end_slab) > size(start_index,2)
    end_slab = num2str(size(start_index,2));                    % update endslab if the data provided is of shorter length
end

n_slices = end_index(str2double(end_slab))-start_index(str2double(start_slab))+1;   % Number of Slices
slices = zeros(pkeyn*skeyn,n_slices,'single');                                      % Initalize matrix for all inlines, xlines and slices

%% Load the data, create slices and Pick the Water Bottom
loopfin = size(job_meta.liveblocks,1);                          % Number of live blocks
lpi = 1;
% check if job meta file has a water bottom path?

% Loop though all live blocks
while lpi <= loopfin
    i_block = job_meta.liveblocks(lpi);                            % Block Number for Current Live Block
    [~, traces, ilxl_read, ~] = ...
        node_segy_read(job_meta_path,vol_index,num2str(i_block));
    %traces = bandpass_filter(traces,job_meta.s_rate/1e6,0,8,100,125);
    traces = [zeros(1,size(traces,2)); traces(2:end,:)];
    % Loop round to pick the water bottom
    %wb_idx = ones(1,size(traces,2));
        wb_idx = water_bottom_picker(traces,10);
        wb_idx(wb_idx < 0) = 1;
        win_sub = bsxfun(@plus,wb_idx,(0:job_meta.n_samples{str2double(vol_index)}-1)');
        win_ind = bsxfun(@plus,win_sub,(0:job_meta.n_samples{str2double(vol_index)}:...
            job_meta.n_samples{str2double(vol_index)}*(size(traces,2)-1)));

      traces(win_sub<=job_meta.n_samples{str2double(vol_index)}) = ...
        traces(win_ind(win_sub<=job_meta.n_samples{str2double(vol_index)}));
    traces(win_sub>job_meta.n_samples{str2double(vol_index)}) = 0;
    [~,n_traces] = size(traces);
    
    % hilbert
    if str2double(hilb) == 1 && str2double(amp) == 0 
        traces = str2double(flip_phase).*imag(hilbert(traces)); % this will perform a -90 phase rotation
    elseif str2double(hilb) == 0 && str2double(amp) == 1 
        traces = sqrt(traces.^2);
    end
    
    % pad traces to account for zeros at end
    pad_start = zeros(floor(ns_win/2),n_traces);
    pad_end = zeros(floor(ns_win/2),n_traces);
    traces = [pad_start;traces;pad_end];
    
    if str2double(start_slab) > size(start_index,2)
        fprintf('Start slab variable too large\n');
        return
    end    
   
    n_slices = end_index(str2double(end_slab))-start_index(str2double(start_slab))+1;
    
    n_iline = (ilxl_read(:,1)-job_meta.pkey_min(str2double(vol_index)))/pkey_inc_mode+1;
    n_xline = (ilxl_read(:,2)-job_meta.skey_min(str2double(vol_index)))/skey_inc_mode+1;
    lin_ind = ((n_iline-1).*skeyn)+n_xline;
    
    if lpi == 1;
        %slices = zeros(pkeyn*skeyn,n_slices,'single');
        max_wb_idx = max(wb_idx);
        max_n_slices = n_slices;
        slices(double(lin_ind),1:n_slices) = traces(start_index(str2double(start_slab)):end_index(str2double(end_slab)),:)';
        win_mid = [start_index;end_index];
    else
        if max_wb_idx < max(wb_idx);
            max_wb_idx = max(wb_idx);
        end
        if n_slices > max_n_slices
            slice_append = zeros(pkeyn*skeyn,n_slices);
            slices = [slices;slice_append];
            slices(double(lin_ind),1:n_slices) = traces(start_index(str2double(start_slab)):end_index(str2double(end_slab)),:)';
            max_n_slices = n_slices;
            win_mid = [start_index;end_index];
        else % write into array
            slices(double(lin_ind),1:n_slices) = traces(start_index(str2double(start_slab)):end_index(str2double(end_slab)),:)';
        end
    end
    fprintf('-- Block %d of %d --\n',lpi,loopfin);
    lpi = lpi + 1;
    
    % Save water bottom pick
    
end
    clear wb_idx_in
    win_mid = [start_index(str2double(start_slab)):start_index(str2double(end_slab));end_index(str2double(start_slab)):end_index(str2double(end_slab))];
    win_mid(3,:) = win_mid(1,:)+floor(ns_win/2);
    n_windows = size(win_mid,2);   
    zero_log = sum(slices,2) == 0;  

%% Calculate anomaly volume
%Loop through all the windows
for ii = 1:n_windows
    block_sz = 100000;                                          % Size of Block
    count = 1;                                                  % Initialize a counter
    data = slices(:,ii:ii+ns_win);                             % Locally load the data for the current window   
    %------------------------------------------------
    % Temp section to store data size
    data = data(:);

    if min(data) ~= max(data)
        data = data(data ~= 0);                                 % remove hard zeros
        data = data(~isnan(data));                              % remove NaN's
        index = (rand(size(data,1),1) <= 0.3);
        data = data(index);
        numdata = length(data);                                 % total number of real datapoints (without zeros and NaN's)
        widthbin = 3.5*std(data)/(numdata^(1/3));               % gives optimum bin width for gaussian distribution (fairly robust for other distributions)
        bins = ((min(data):widthbin:max(data)))';               % defines bin centres based on bin width and data extremes
        data = interp1(bins,bins,data,'nearest','extrap');      % interpolates data onto closest bin center
        check_bins = ~ismember(bins,data);                      % finds empty bins
        data = sort([data;bins(check_bins)]);                   % temporarily pads data with empty bins to ensure no bins are empty, and sorts low to high
        [~,datapdf,~] = unique(data);                           % returns the highest index of the unique values in data
        datapdf = diff([0;datapdf]);                            % difference between indices is the number of values in a bin
        datapdf(check_bins) = 0;                                % reset the padded bins to zero
        %datapdf = datapdf_tmp;                                 % put into a cell array
        datacdf = cumsum(datapdf)-datapdf/2;                    % cumulative sum across bins and re-centers to bin center
        datacdf = datacdf/max(datacdf);                         % normalises max of CDF to 1
        
        % Calculate probabilities from cdf
        % get central slice in window
%         if str2double(amp) == 1
%             data = interp1(bins,bins,sqrt(slices(:,ii+floor(ns_win/2)).^2),'nearest','extrap');         % interpolates data onto closest bin center
%             [~,idxdataprob] = ismember(data,bins);                                                      % finds the bin index that the data is in        
%             dataprob = datacdf(idxdataprob);
%         else
            data = interp1(bins,bins,slices(:,ii+floor(ns_win/2)),'nearest','extrap');                  % interpolates data onto closest bin center
            [~,idxdataprob] = ismember(data,bins);                                                      % finds the bin index that the data is in        
            dataprob = 1-datacdf(idxdataprob);
        %end

        dataprob(zero_log,1) = 0;
    else
        fprintf('Blank slice at %d\n',ii);
        dataprob = zeros(pkeyn*skeyn,1);
    end
        
    for ii_block = 1:block_sz:size(dataprob,1)
        maxblock = ii_block + block_sz - 1;
        
        if maxblock >  size(dataprob,1)
            maxblock = size(dataprob,1);
            block_sz = (maxblock - ii_block+1);
        end
        
        fid = fopen([job_meta.output_dir,'p_',num2str(count),'_sas_result_slice_il',num2str(pkeyn),'_xl',num2str(skeyn),'_',num2str(win_mid(1,1)),'-',...
            num2str(win_mid(1,end)),'_traces',num2str(block_sz),'_block.bin'],'a');
        
        fwrite(fid,dataprob(ii_block:maxblock,1),'float32');
        fclose(fid);
        count = count + 1;
    end
    fprintf('Completed %d of %d\n',ii,n_windows);
end

%fwrite(fid,dataprob','float32');

end
