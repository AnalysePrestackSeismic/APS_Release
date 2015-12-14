
% residual_velocity_analysis('/URY/segy/2014_BG_water_column_imaging/matlabjob_meta/job_meta_23Jul2014.mat','10','velocity_path');

% 'velocity_path' parameter is not used and is only there for future
% implementation

total_picks=0;
total_traces=0;

for iblock = 2:23
    matfile = strcat('/segy/URY/2014_BG_water_column_imaging/matlab/rms_picks_test_block',int2str(iblock),'.mat');
    load(matfile,'ntraces','tvpairs');
    total_traces = total_traces + ntraces;
    
    for trace = 1:ntraces
        % copy the first pick to the first time to avoid NaNs in gridding
        % copy last pick to end for interval vel conversion
        current_trace = trace+total_traces-ntraces;
        numpicks(current_trace)=size(tvpairs{trace},1)+2; % add 2 to numpicks for first/last padding
        pick_idx(current_trace)=total_picks+1;
        total_picks=total_picks+numpicks(current_trace);
        
        alltimes(pick_idx(current_trace))=1;
        allvels(pick_idx(current_trace))=tvpairs{trace}(1,2);
        
        alltimes(total_picks)=5000;
        allvels(total_picks)=tvpairs{trace}(end,2);
        
        alltimes(pick_idx(current_trace)+1:total_picks-1)=tvpairs{trace}(:,1);
        allvels(pick_idx(current_trace)+1:total_picks-1)=tvpairs{trace}(:,2);
        alllocs(pick_idx(current_trace):total_picks)=current_trace;
    end
    
end


% figure;scatter(alllocs,-alltimes,50,allvels,'filled'); caxis([1.48 1.54]);

% triangulate onto a grid
[binx,biny]=meshgrid([1:total_traces],[1:5000]);
gridtv=griddata(alllocs,alltimes,allvels,binx,biny,'linear');

save('/segy/URY/2014_BG_water_column_imaging/matlab/grid_rms_picks.mat');


        