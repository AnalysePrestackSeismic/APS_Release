
load('/segy/URY/2014_BG_water_column_imaging/matlab/grid_rms_picks.mat');

% test whether smoothing in slowness makes any difference

gridtv=gridtv.^-1;

% gaussian smoothing
smth_grid_tv=gaussian_1dsmth(gridtv,4001);

% vertical smoothing

smth_grid_tv=(gaussian_1dsmth(smth_grid_tv',201))';

smth_grid_tv=smth_grid_tv.^-1;

% extract smoothed vels at picked locations
for ii=1:total_picks
    allsmvels(ii)=smth_grid_tv(round(alltimes(ii)),alllocs(ii));
end

% figure; scatter(alllocs,-alltimes,20,allvels,'filled'); caxis([1.48 1.54]);
% figure; scatter(alllocs,-alltimes,20,allsmvels,'filled'); caxis([1.48 1.54]);

% interval velocities

allvint=zeros(size(allvels));

for trace = 1:total_traces;
    
    ttimes = alltimes(pick_idx(trace):pick_idx(trace)+numpicks(trace)-1);
    vrms = allsmvels(pick_idx(trace):pick_idx(trace)+numpicks(trace)-1);
    
    tts = [0 ttimes(1:end-1)];
    vrms_s = [vrms(1) vrms(1:end-1)];
    
    allvint(pick_idx(trace):pick_idx(trace)+numpicks(trace)-1) = sqrt(((vrms.^2.*ttimes) - (vrms_s.^2.*tts)) ./ (ttimes-tts));
    
end

figure; scatter(alllocs,-alltimes,20,allvint,'filled'); caxis([1.48 1.54]);

% create interval vel traces for display

% copy each pick to immediately after previous pick
% then interp1d to fill between
% and reshape to make into a grid

alltimes2 = zeros(size(alltimes,2)*2-1,1);
allvels2 = zeros(size(allvint,2)*2-1,1);

alltimes2(1:2:end)=round(alltimes+(alllocs-1).*5100);
alltimes2(2:2:end-1)=round(alltimes(1:end-1)+1+(alllocs(1:end-1)-1).*5100);

allvels2(1:2:end)=allvint;
allvels2(2:2:end-1)=allvint(2:end);


grid_vint = reshape(interp1(alltimes2,allvels2,1:5100*total_traces),5100,total_traces);

grid_vint10 = gaussian_1dsmth(grid_vint,21);

grid_vint10 = grid_vint10(:,10:10:end);

