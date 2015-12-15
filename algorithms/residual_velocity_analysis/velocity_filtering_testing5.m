
load('/data/URY/segy/2014_BG_water_column_imaging/conditioned_datasets/cdps/4150A216/gridfit_rms_picks_800x40.mat');

% for very shallow, calculate interval vels sample by sample because don't
% have picks

% merge that with interval vels at the pick locations for rest of section


% interval vels at pick locations =========================================

% decimate original pick times onto new grid

alllocs10=alllocs(mod(alllocs,10)==mod(dec_ilxl(1,2),10));
alltimes10=alltimes(mod(alllocs,10)==mod(dec_ilxl(1,2),10));
numpicks10=size(alllocs10,2);

% extract smoothed vels at picked locations
alllocs_grid = ((alllocs10-dec_ilxl(1,2))./10)+1;

for ii=1:numpicks10
    allsmvels10(ii)=gridtv(ceil(alltimes10(ii)/10),alllocs_grid(ii));

end

% figure; scatter(alllocs,-alltimes,20,allvels,'filled'); caxis([1.48 1.54]);
% figure; scatter(alllocs,-alltimes,20,allsmvels,'filled'); caxis([1.48 1.54]);

% interval velocities

allvint10=zeros(size(allsmvels10));

for trace = dec_ilxl(:,2)';
    picks=alllocs10==trace;
    ttimes=alltimes10(picks);
    vrms=allsmvels10(picks);
    tts = [0 ttimes(1:end-1)];
    vrms_s = [vrms(1) vrms(1:end-1)];
    allvint10(picks) = sqrt(((vrms.^2.*ttimes) - (vrms_s.^2.*tts)) ./ (ttimes-tts));
    
end

% figure; scatter(alllocs,-alltimes,20,allvint,'filled'); caxis([1.48 1.54]);

% create interval vel traces for display

% copy each pick to immediately after previous pick
% then interp1d to fill between
% and reshape to make into a grid


alltimes10_2 = zeros(size(alltimes10,2)*2-1,1);
allvels10_2 = zeros(size(allvint10,2)*2-1,1);

alltimes10_2(1:2:end)=round(alltimes10+(alllocs_grid-1).*5100);
alltimes10_2(2:2:end-1)=round(alltimes10(1:end-1)+1+(alllocs_grid(1:end-1)-1).*5100);

allvels10_2(1:2:end)=allvint10;
allvels10_2(2:2:end-1)=allvint10(2:end);

dec_traces = size(dec_ilxl,1);

grid_vint = reshape(interp1(alltimes10_2,allvels10_2,1:5100*dec_traces),5100,dec_traces);

% gridded dix interval velocities for shallow =============================

vint_shallow=zeros(50,dec_traces);

for trace = 1:dec_traces;
    
    ttimes = (10:10:500)';
    vrms = gridtv(1:50,trace);
    
    tts = [0;ttimes(1:end-1)];
    vrms_s = [vrms(1);vrms(1:end-1)];
    
    vint_shallow(:,trace) = sqrt(((vrms.^2.*ttimes) - (vrms_s.^2.*tts)) ./ (ttimes-tts));
    
    vint_merge(:,trace) = interp1([10:10:200 250:5000],[vint_shallow(1:20,trace)' grid_vint(250:5000,trace)'],1:5000,'linear','extrap');
end






