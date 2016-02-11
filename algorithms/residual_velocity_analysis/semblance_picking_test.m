
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
semb_meta_path = '/data/URY/segy/2014_BG_water_column_imaging/conditioned_datasets/semblance/4110A212/job_meta/job_meta_06Nov2014.mat';

semb_meta = load(semb_meta_path);
ntraces = ((semb_meta.skey_max-semb_meta.skey_min)/semb_meta.skey_inc)+1;
peak_semb = zeros(semb_meta.n_samples{1},ntraces);
peak_semb_thr = zeros(semb_meta.n_samples{1},ntraces);
pick_vel = zeros(semb_meta.n_samples{1},ntraces);


for block = semb_meta.liveblocks'
    
    block_idx = find(semb_meta.liveblocks==block);    
    [semblance, traces, il_xl, offsets] = node_segy_read(semb_meta_path,'1',num2str(block));
    fxl = semb_meta.block_keys(block_idx,3);
    lxl = semb_meta.block_keys(block_idx,4);
    fcdp = ((fxl-semb_meta.skey_min)/semb_meta.skey_inc)+1;
    lcdp = ((lxl-semb_meta.skey_min)/semb_meta.skey_inc)+1;
   
    semblance.fold = max(semblance.trace_ilxl_bytes(:,7));
    
    traces = reshape(traces,size(traces,1),semblance.fold,[]);
    
    semblance.n_gathers = size(traces,3);
    
    threshold = 12000; gap = 25; med_length=5; med_threshold = 10;
       
    
    for cdp_num = fcdp:lcdp % cdp_num is the index of the cdp along the whole line
        gather_idx = cdp_num-fcdp+1; % gather_idx is the index within the block
        semb = traces(:,:,gather_idx);
        [peak_semb(:,cdp_num),vel_idx] = max(semb,[],2); % pick peak semblance at every time
        peak_semb_copy = peak_semb(:,cdp_num);
        peak_semb_copy(1:100)=0; % zero first 100 samples
        
        while max(peak_semb_copy)>threshold
            
            [peak_val,peak_idx] = max(peak_semb_copy);
            peak_semb_thr(peak_idx,cdp_num) = peak_val;
            pick_vel(peak_idx,cdp_num) = vel_idx(peak_idx);
            
            start_mask=max(peak_idx-gap,1);
            end_mask=min(peak_idx+gap,semblance.n_samples);
            
            peak_semb_copy(start_mask:end_mask) = 0;
            
        end
        
    end
    
end


[semb_times,semb_locs]=find(peak_semb_thr);
semb_vels = pick_vel(peak_semb_thr>0).*2+1448;

% figure; scatter(semb_locs,semb_times.*-1,20,semb_vels,'filled'); caxis([1490 1540]);

    % median filter to get rid of anomalous values
    % ignoring the fact that the vels aren't regularly sampled


 for trace = 1:ntraces

     vels = semb_vels(semb_locs==trace);
     times = semb_times(semb_locs==trace);
     vels_filt = medfilt1(vels,med_length);
     vels_error = vels - vels_filt;
     vels(vels_error>med_threshold) = vels_filt(vels_error>med_threshold);
     semb_vels(semb_locs==trace) = vels;
 end
 
figure; scatter(semb_locs,semb_times.*-1,20,semb_vels,'filled'); caxis([1490 1540]);

gridvel=gridfit(semb_locs,semb_times,semb_vels,1:19315,10:10:2500,'smoothness',[1600 200]);

