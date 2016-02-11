function [] = segy_plot_ava(job_meta_path,s_il,e_il,s_xl,e_xl,il_inc,xl_inc,startvol,volinc,endvol,trace_decimate,start_time,end_time)
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
%maxzout = str2double(maxzout);

% 1KG
% 

job_meta = load(job_meta_path);
live_block_keys = [job_meta.block_keys(job_meta.liveblocks,:) job_meta.liveblocks];

s_il = strsplit(s_il,' ');
e_il = strsplit(e_il,' ');
s_xl = strsplit(s_xl,' ');
e_xl = strsplit(e_xl,' ');

n_pos = size(s_il,2);
for i_pos = 1:1:n_pos
    s_il{i_pos} = str2num(s_il{i_pos});
    e_il{i_pos} = str2num(e_il{i_pos});
    s_xl{i_pos} = str2num(s_xl{i_pos});
    e_xl{i_pos} = str2num(e_xl{i_pos});
end

s_il = cell2mat(s_il);
e_il = cell2mat(e_il);
s_xl = cell2mat(s_xl);
e_xl = cell2mat(e_xl);

for i_pos = 1:1:n_pos
    log_block = live_block_keys(:,1) <= s_il(i_pos) & live_block_keys(:,2) >= e_il(i_pos) & live_block_keys(:,3) <= s_xl(i_pos) & live_block_keys(:,4) >= e_xl(i_pos);
    block_keep{:,i_pos} = live_block_keys(log_block,5);
end
startvol = str2double(startvol);    % number of the first angle trace/volume to read
volinc = str2double(volinc);        % angle trace/volume increment 
%endvol = job_meta.nvols;
endvol = str2double(endvol);        % number of the last angle trace/volume to read
totalvol = length(startvol:volinc:endvol);  

il_inc = str2double(il_inc); 
xl_inc = str2double(xl_inc); 

padding = 50;       % water bottom pick padding
extrapad = 0;     % extra padding to top of dataset, this many samples will get zeroed below the wb 
start_sample = floor(str2double(start_time)/(job_meta.s_rate/1000));
end_sample = floor(str2double(end_time)/(job_meta.s_rate/1000));
maxzout = end_sample;

%%
for i_pos = 1:1:n_pos
    for i_block = 1:1:size(block_keep{:,i_pos})
        [~, vol_traces, ilxl_read, offset_read] = node_segy_read(job_meta_path,'1',num2str(block_keep{i_block,i_pos}));  
        %uniq_ilxl_read{i_pos,i_block} = unique(ilxl_read{i_pos,i_block},'rows');
        
        log_pos = ilxl_read(:,1) >= s_il(i_pos) & ilxl_read(:,1) <= e_il(i_pos) & ilxl_read(:,2) >= s_xl(i_pos) & ilxl_read(:,2) <= e_xl(i_pos);

        traces{i_pos,i_block} = vol_traces(start_sample:end_sample,log_pos);
        ilxl_pos{i_pos,i_block} =  ilxl_read(log_pos,:);
        offset_pos{i_pos,i_block} = offset_read(log_pos);
    end 
    
end

%% Plot gathers


subplot(2,1,1); imagesc(traces{1},[-650000 650000]); subplot(2,1,2); imagesc(traces{2},[-650000 650000]);
colormap('gray')

xl_1 = 9084;
xl_2 = 9104;

amp{1} = traces{1,1}(:,ilxl_pos{1,1}(:,2) == xl_1);
input_angles{1} = offset_pos{1,1}(:,ilxl_pos{1,1}(:,2) == xl_1);

amp{2} = traces{2,1}(:,ilxl_pos{2,1}(:,2) == xl_2);
input_angles{2} = offset_pos{2,1}(:,ilxl_pos{1,1}(:,2) == xl_2);

x_d = [input_angles{1}(1,1):1:input_angles{1}(1,end)];
y_d = [str2double(start_time):5:str2double(end_time)];
plot_1_s = 241;
plot_1_e = 251;

plot_2_s = 477;
plot_2_e = 487;

mid_p1 = round(plot_1_s+(plot_1_e-plot_1_s)/2)*5+3000;
mid_p2 = round(plot_2_s+(plot_2_e-plot_2_s)/2)*5+3000;

figure(2)
subplot(3,2,1); imagesc(x_d,y_d,amp{1}); title(sprintf('Mzia il: %d xl: %d \n',s_il(1),xl_1)); hold all; plot(input_angles{1},repmat(mid_p1,1,size(input_angles{1},2))); hold off; ylabel('Depth (m)'); xlabel('Angle (Degrees)'); 
subplot(3,2,3); plot((sind(double(input_angles{1}))).^2,amp{1}(plot_1_s:plot_1_e,:),...    '--rs',...
                            'LineWidth',2,...
                            'MarkerEdgeColor','k',...
                            'MarkerFaceColor','g',...
                            'MarkerSize',10)
                        xlabel('sin^2(angle)');
subplot(3,2,2); imagesc(x_d,y_d,amp{2}); title(sprintf('1BG il: %d xl: %d \n',s_il(2),xl_2)); hold all; plot(input_angles{1},repmat(mid_p2,1,size(input_angles{2},2))); hold off; ylabel('Depth (m)'); xlabel('Angle (Degrees)'); 
subplot(3,2,4); plot((sind(double(input_angles{2}))).^2,amp{2}(plot_2_s:plot_2_e,:),...    '--rs',...
                            'LineWidth',2,...
                            'MarkerEdgeColor','k',...
                            'MarkerFaceColor','g',...
                            'MarkerSize',10)
                        xlabel('sin^2(angle)');  
y_zoom1 = [(plot_1_s+1)*5+3000:5:(plot_1_e+1)*5+3000];
y_zoom2 = [(plot_2_s+1)*5+3000:5:(plot_2_e+1)*5+3000];
subplot(3,2,5); imagesc(x_d,y_zoom1,amp{1}(((y_zoom1-3000)/5)+1,:)); ylabel('Depth (m)'); xlabel('Angle (Degrees)');
subplot(3,2,6); imagesc(x_d,y_zoom2,amp{2}(((y_zoom2-3000)/5)+1,:)); ylabel('Depth (m)'); xlabel('Angle (Degrees)');
colormap('gray');

imagesc
imagesc(traces{1,1}(:,ilxl_pos{1,1}(:,2) == 6340))
    figure(2)
    subplot(ceil(sqrt(n_samples)),ceil(sqrt(n_samples)),i_sample);
    for i_trace=1:str2double(trace_decimate):n_traces
        if i_trace == 1
            plot((sind(input_angles)).^2,gather(:,i_sample,i_trace),...
                            '--rs',...
                            'LineWidth',2,...
                            'MarkerEdgeColor','k',...
                            'MarkerFaceColor','g',...
                            'MarkerSize',10)
        else
            hold on
            plot((sind(input_angles)).^2,gather(:,i_sample,i_trace),...
                '--rs',...
                'LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10)
            hold off
        end
    end
    title(['Sample ',num2str(i_sample)])            
    xlabel('sin^2(\Theta)')
    ylabel('Amplitude')  
    if i_sample == n_samples
        subplot(ceil(sqrt(n_samples)),ceil(sqrt(n_samples)),i_sample+1);
        imagesc(gather(:,:,i_trace)');
    end


vol_index_wb = 1;
if job_meta.is_gather == 0
    pick_wb_ind = ceil(job_meta.nvols*0.6667);
    
    [~, traces{vol_index_wb}, ilxl_read{vol_index_wb}] = node_segy_read(job_meta_path,num2str(pick_wb_ind),i_block);
    % check to make sure it read something if not exit
    if size(traces{vol_index_wb},1) == 1 && size(traces{vol_index_wb},2) == 1
        return
    end
    traces{vol_index_wb} = traces{vol_index_wb}(1:maxzout,:);
    n_samples_fullz = job_meta.n_samples{1};
    job_meta.n_samples{1} = maxzout;    
else
    
    
    % find the water bottom on a few gathers and take the middle one t use
    % as an offset plane to pick the water bottom on
    
    % read all the data for this block
    % node_segy_read(job_meta_path,vol_index,i_block)
    [~, vol_traces, ilxl_read{vol_index_wb}, offset_read] = node_segy_read(job_meta_path,'1',i_block);    
    vol_traces = vol_traces(1:maxzout,:);       % truncate data to make z axis number
    job_meta.n_samples{1} = maxzout;            % Write max zout in job meta file    
    offset = unique(offset_read);               % find the total number of offsets
    %#############for lobster TRN limit offset to 40 and rerun for failed
    %blocks###################################################################
    
    if isempty(vol_traces) == 1 &&  isempty(ilxl_read) == 1 && isempty(offset_read) == 1
        return
    end
    %juststack = 1;
    %if juststack == 1
    traces{vol_index_wb} = zeros(size(vol_traces,1),size(vol_traces,2)/size(offset,2));
    for stki  = 1:size(offset,2)
        traces{vol_index_wb} = traces{vol_index_wb} + vol_traces(:,offset_read == offset(stki));
    end
    
    pick_wb_ind = job_meta.live_offset_avg;        
end

% Pick water bottom or use a pre picked water bottom horizon
if isfield(job_meta, 'wb_path')
    %wb_idx = dlmread(job_meta.wb_path,'delimiter','\t');
    wb_idx_in = dlmread(job_meta.wb_path);
    %wb_idx_in = sortrows(wb_idx_in,[1 2]);
    % col 1 inline
    % col 2 xline
    % col 3 twt
    if job_meta.is_gather == 1
        [~,locations] = ismember(ilxl_read{1}(1:length(offset):end,:),wb_idx_in(:,1:2),'rows');
    else    
        [~,locations] = ismember(ilxl_read{1}(1:end,:),wb_idx_in(:,1:2),'rows');    
    end
    %wb_idx = zeros(size(traces{vol_index_wb},2),1);
    zero_loc = locations ~= 0;
    %wb_idx
    %zero_loc = wb_idx_in(locations(zero_loc),3); comment this line to
    %match wavelet estimation - cj 27_02_2015
    % but i do not see why this is not just
    %wb_idx =  wb_idx_in(locations,3)
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
    
    wb_idx = round(wb_idx-padding);
    wb_idx(isnan(wb_idx)) = 1;
    wb_idx(wb_idx < 1) = 1;
    win_sub = bsxfun(@plus,wb_idx,(0:job_meta.n_samples{vol_index_wb}-max(wb_idx))');
    
    win_ind = bsxfun(@plus,win_sub,(0:job_meta.n_samples{vol_index_wb}:...
        job_meta.n_samples{vol_index_wb}*(size(traces{vol_index_wb},2)-1)));
else
    [wb_idx] = water_bottom_picker(traces{vol_index_wb},padding);
    wb_idx(wb_idx < 0) = 1;
    [max_wb,max_wb_ind] = max(wb_idx);
    wb_idx(max_wb_ind) = 0;
    [max_wb2,max_wb_ind2] = max(wb_idx);
    
    if max_wb - max_wb2 > 20
        max_wb_idx = max_wb2;
        wb_idx(max_wb_ind) = max_wb2;
    else
        wb_idx(max_wb_ind) = max_wb;
        max_wb_idx = max_wb;
    end
    
    win_sub = bsxfun(@plus,wb_idx,(0:job_meta.n_samples{vol_index_wb}-max_wb_idx)');   
    
    
    win_ind = bsxfun(@plus,win_sub,(0:job_meta.n_samples{vol_index_wb}:...
    job_meta.n_samples{vol_index_wb}*(size(traces{vol_index_wb},2)-1)));
   
end

% read in all the data and decimate as required (allready read in for gathers)
if job_meta.is_gather == 1
    % reshape the gather array to the same 3d matrix as the angle volumes
    %vol_tracescjtmp = reshape(vol_traces,size(traces{vol_index_wb},1),length(offset),size(traces{vol_index_wb},2));
    %vol_traces = vol_tracescjtmp(:,(startvol:volinc:endvol),:);
    %clear vol_tracescjtmp;
    % this can just be vol_traces_cjtmp is not needed, just put in to check
    % debug stages
    vol_traces = reshape(vol_traces,size(traces{vol_index_wb},1),length(offset),size(traces{vol_index_wb},2));
    
    % try applying a running average to the data to smooth out variations
    % on the angle gathers
        
    filttraces = [1 1 1]/3;
    %----- Finding the mute and discarding the points outside it by a mask----------------
    for ii = 1:1:size(vol_traces,3)
        % mute the low amplitude areas
        mask = low_amp_mute(vol_traces(:,:,ii));
        % apply a small smoother to reduce noise
%        vol_traces(:,:,ii) = conv2(1,filttraces,vol_traces(:,:,ii),'same');
        vol_traces(:,:,ii) = vol_traces(:,:,ii) .* mask;
        %imagesc(vol_traces(:,:,ii));
        %colormap(gray);
        %caxis([-500 500]);
    end
    
    %drop the angles that are not needed
    % note the 3 trace average above to allow 2:1 decimation with a bit of
    % anti aliasing, not perfect but ...
    vol_traces = vol_traces(:,(startvol:volinc:endvol),:);    
    
    input_angles = double(offset(startvol:volinc:endvol));
    % resize the ilxl data to drop to stack fold
    tmp_ilxlrd = ilxl_read{vol_index_wb}(1:length(offset):end,:);
    ilxl_read{vol_index_wb} = tmp_ilxlrd;
    clear tmp_ilxlrd;
else
    % Load block for remaining angle traces
    % initialise the data array
    vol_traces = zeros(totalvol,n_samples_fullz,size(traces{vol_index_wb},2));
    
    % if tottracerun == 0;
    %     vol_traces = zeros(totalvol,size(traces{vol_index_wb},1),size(traces{vol_index_wb},2));
    % else
    %     vol_traces = zeros(totalvol,size(traces{vol_index_wb},1),tottracerun);
    % end

    vol_count = 1;
    for i_vol = startvol:volinc:endvol
        fprintf('reading data no of vol = %d\n',i_vol)

        % Read traces
        %[~, traces{vol_count}, ~, ~] = node_segy_read(job_meta_path,num2str(i_vol),i_block);
        [~, vol_traces(vol_count,:,:), ~, ~] = node_segy_read(job_meta_path,num2str(i_vol),i_block);

        % Flatten traces to water bottom
        %traces{vol_count} = traces{vol_count}(win_ind);
        %vol_traces(vol_count,:,:) = vol_traces(vol_count,win_ind);
        % truncate the dataset if just running a test for a set number of
        % traces
        %     if tottracerun ~= 0;
        % %         if vol_count == 1;
        % %             vol_traces = zeros(totalvol,size(traces{vol_count},1),tottracerun);
        % %         end
        %         %###########################################################
        %         % cj test edit resize the dataset temp to test with
        %         vol_traces(vol_count,:,:) = traces{vol_count}(:,1:tottracerun);
        %         traces{vol_count} = traces{vol_count}(:,1:tottracerun);
        %     end
        
        % keep a list of the angle values read in, this might need to change
        % for angle gathers with on a single angle value; it needs to be
        % assigned to both fields in the job meta file
        input_angles(vol_count) = (job_meta.angle{i_vol}(2)+job_meta.angle{i_vol}(1))/2;
        
        if plot_on == 2
            figure(1)
            subplot(1,totalvol,vol_count); imagesc(vol_traces(vol_count,:,:));
        end
        
        vol_count = vol_count + 1;
    end

    % deal with zeros in the IGmatrix build rather than as a mask as
    % already zeros, this is just to cover low amp noise
    %----- Finding the mute and discarding the points outside it by a mask----------------
    for ii = 1:1:size(vol_traces,3)
        % mute the low amplitude areas
        mask = low_amp_mute(vol_traces(:,:,ii)');
        % apply a small smoother to reduce noise
        %vol_traces(:,:,ii) = conv2(1,filttraces,vol_traces(:,:,ii),'same');
        vol_traces(:,:,ii) = vol_traces(:,:,ii) .* mask';
        
        %imagesc(vol_traces(:,:,ii));
        %colormap(gray);
        %caxis([-500 500]);
    end    

    % truncate data to make z axis number
    vol_traces = vol_traces(:,1:maxzout,:);
    job_meta.n_samples{1} = maxzout;
end
clear traces win_sub wb_idx;
%%
%==========================================================================================================================
% now flatten to the water bottom offset plane at a time
%
padding = padding + extrapad;
%
if job_meta.is_gather == 1
    for i_vol = 1:totalvol
        tmp_vol = squeeze(vol_traces(:,i_vol,:));
        vol_traces(1:size(win_ind,1),i_vol,:) = tmp_vol(win_ind);
    end
    % resize the traces to drop the trailing blanks after the shift to the wb
    vol_traces = vol_traces(1:size(win_ind,1),:,:);
    % taper the top 25 samples to avoid wb reflection creating artifacts due to
    % high amplitude

    % need to apply a taper that does not go to zero, make the wavelet try
    % to match zero and does not work
    %taper = single([ zeros(floor(padding-5),1)',linspace(0,1,20)])';
    %vol_traces(1:length(taper),:,:) = bsxfun(@times,vol_traces(1:length(taper),:,:),taper);
    
else
    for i_vol = 1:totalvol
        tmp_vol = squeeze(vol_traces(i_vol,:,:));
        % flattern to water bottom
        vol_traces(i_vol,1:size(win_ind,1),:) = tmp_vol(win_ind);
    end
    % resize the traces to drop the trailing blanks after the shift to the wb
    vol_traces = vol_traces(:,1:size(win_ind,1),:);
    % taper the top 25 samples to avoid wb reflection creating artifacts due to
    % high amplitude
    
    % need to apply a taper that does not go to zero, make the wavelet try
    % to match zero and does not work
    %taper = single([ zeros(floor(padding-5),1)',linspace(0,1,20)]);
    %vol_traces(:,1:length(taper),:) = bsxfun(@times,vol_traces(:,1:length(taper),:),taper);
    
end

for i_vol = 1:1:job_meta.nvols 
    if i_vol ~= str2double(vol_index_wb) % don't repeat load for previously read stack
        % Read traces
        [~, traces{i_vol}, ~, ~] = ...
        node_segy_read(job_meta_path,num2str(i_vol),i_block);
    end   
    % Flatten traces to water bottom
    traces{i_vol} = traces{i_vol}(win_ind);  
    gather(i_vol,1:n_samples,1:n_traces) = traces{i_vol}(start_sample:end_sample,:);
    %traces{i_vol} = traces{i_vol}(,:);
    input_angles(i_vol) = (job_meta.angle{i_vol}(2)+job_meta.angle{i_vol}(1))/2;

    figure(1)
    subplot(1,job_meta.nvols,i_vol); imagesc(traces{i_vol});   
    hold all
    plot(repmat(start_sample,1,n_traces))
    plot(repmat(end_sample,1,n_traces))
    hold off
end

% Load DIGI result

for i_sample = 1:1:n_samples
    figure(2)
    subplot(ceil(sqrt(n_samples)),ceil(sqrt(n_samples)),i_sample);
    for i_trace=1:str2double(trace_decimate):n_traces
        if i_trace == 1
            plot((sind(input_angles)).^2,gather(:,i_sample,i_trace),...
                            '--rs',...
                            'LineWidth',2,...
                            'MarkerEdgeColor','k',...
                            'MarkerFaceColor','g',...
                            'MarkerSize',10)
        else
            hold on
            plot((sind(input_angles)).^2,gather(:,i_sample,i_trace),...
                '--rs',...
                'LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10)
            hold off
        end
    end
    title(['Sample ',num2str(i_sample)])            
    xlabel('sin^2(\Theta)')
    ylabel('Amplitude')  
    if i_sample == n_samples
        subplot(ceil(sqrt(n_samples)),ceil(sqrt(n_samples)),i_sample+1);
        imagesc(gather(:,:,i_trace)');
    end
end
% Test unflatten
% Unflatten data
% [ns,ntraces] = size(traces{1});
% traces_unflat = [traces{1};zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
% for kk = 1:length(wb_idx)
%     traces_unflat(:,kk) = circshift(traces_unflat(:,kk),wb_idx(kk));
% end
%% Plot AVA