function [semb_times semb_locs semb_vels semb_filt_vels] = pick_semblance(semb_meta_path)

% semb_meta_path = '/data/URY/segy/2014_BG_water_column_imaging/conditioned_datasets/semblance/4110A212/job_meta/job_meta_06Nov2014.mat';

semb_meta = load(semb_meta_path);

ntraces = ((semb_meta.skey_max-semb_meta.skey_min)/semb_meta.skey_inc)+1;
peak_semb = zeros(semb_meta.n_samples{1},ntraces);
peak_semb_thr = zeros(semb_meta.n_samples{1},ntraces);
pick_vel = zeros(semb_meta.n_samples{1},ntraces);
all_xls = zeros(ntraces,1);

for block = semb_meta.liveblocks'
    block_idx = find(semb_meta.liveblocks==block);
    [semblance, traces, il_xl, offsets] = node_segy_read(semb_meta_path,'1',num2str(block));
    block_xls = unique(il_xl(:,2));
    fxl = block_xls(1); lxl = block_xls(end);
    fcdp = ((fxl-semb_meta.skey_min)/semb_meta.skey_inc)+1;
    lcdp = ((lxl-semb_meta.skey_min)/semb_meta.skey_inc)+1;
    
    semblance.fold = max(semblance.trace_ilxl_bytes(:,7));
    
    traces = reshape(traces,size(traces,1),semblance.fold,[]);
    
    semblance.n_gathers = size(traces,3);
    
    threshold = 10000; gap = 25; med_length=5; med_threshold = 10;
    
    
    for cdp_num = 1:size(block_xls,1) % loop through the cdps in the block
        xl = block_xls(cdp_num); % xl is the actual xl value
        skey = ((xl - semb_meta.skey_min)/semb_meta.skey_inc)+1; % skey is the index from first to last xl on line
        all_xls(skey) = xl;
        semb = traces(:,:,cdp_num);
        [peak_semb(:,skey),vel_idx] = max(semb,[],2); % pick peak semblance at every time
        peak_semb_copy = peak_semb(:,skey);
        peak_semb_copy(1:100)=0; % zero first 100 samples
        
        % keep max semblance picks above threshold with minimum gap between
        % them
        
        while max(peak_semb_copy)>threshold
            
            [peak_val,peak_idx] = max(peak_semb_copy);
            peak_semb_thr(peak_idx,skey) = peak_val;
            pick_vel(peak_idx,skey) = vel_idx(peak_idx);
            
            start_mask=max(peak_idx-gap,1);
            end_mask=min(peak_idx+gap,semblance.n_samples);
            
            peak_semb_copy(start_mask:end_mask) = 0;
            
        end
        
    end
    
end


[semb_times,semb_locs]=find(peak_semb_thr);
semb_locs = all_xls(semb_locs);
% semb_locs = ((semb_locs-1).*semb_meta.skey_inc)+semb_meta.skey_min; % convert index to xl number
semb_vels = pick_vel(peak_semb_thr>0).*2+1448;
semb_filt_vels = zeros(size(semb_vels));

% median filter to get rid of anomalous values
% ignoring the fact that the vels aren't regularly sampled


semb_filt_vels = semb_vels; % ***** license for medfilt1 not available

% for xl = all_xls';
%     pick_idx = semb_locs==xl;
%     vels = semb_vels(pick_idx);
%     times = semb_times(pick_idx);
%     vels_filt = medfilt1(vels,med_length);
%     vels_error = vels - vels_filt;
%     vels(vels_error>med_threshold) = vels_filt(vels_error>med_threshold);
%     semb_filt_vels(pick_idx) = vels;
% end
% 
% figure; scatter(semb_locs,semb_times.*-1,20,semb_vels,'filled'); caxis([1490 1540]);

% take the gridding out of this function because it needs more testing

% fgridxl = round(semb_meta.skey_min/10)*10;
% lgridxl = round(semb_meta.skey_max/10)*10;
% 
% rmsgrid=gridfit(semb_locs,semb_times,semb_vels,fgridxl:10:lgridxl,10:10:2500,'smoothness',[1600 200]);

end
