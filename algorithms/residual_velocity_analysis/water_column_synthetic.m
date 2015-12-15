

vsp9 = importdata('/data/URY/segy/2014_BG_water_column_imaging/TS_Dips/vsp9_descent.txt');

[wavelet, wavelet_times] = ricker(0.002,35);


% slowness = vsp9.data(:,1).^-1;
% depth2 = [0;vsp9.data(1:end-1,2)];  
% ttimes = 2*cumsum(slowness.*(vsp9.data(:,2)-depth2),1);

%depth_time = [vsp9.data(:,2) vsp9.data(:,2)];

[~,sort_index] = sort(vsp9.data(:,2));

vsp_sort_depths = vsp9.data(sort_index,2);
vsp_sort_density = vsp9.data(sort_index,6);
vsp_sort_velocity = vsp9.data(sort_index,1);

[~,sort_index_uniq,~] = unique(vsp_sort_depths);

vsp_depths = vsp_sort_depths(sort_index_uniq);
vsp_density = vsp_sort_density(sort_index_uniq);
vsp_velocity = vsp_sort_velocity(sort_index_uniq);

depths = [2:2:3200]';


density = interp1(vsp_depths,vsp_density,depths);
velocity = interp1(vsp_depths,vsp_velocity,depths);

slowness = velocity.^-1;
depths2 = [0;depths(1:end-1)];  
ttimes = 2*cumsum(slowness.*(depths-depths2),1);

depth_time = [depths ttimes];

% create synthetic traces
[synthetic,syntimes,rc,~,~] = seismogram(velocity,density,depths,wavelet,wavelet_times,depth_time,'',0,1);

% make vrms velocity
[vsp_times ~] = vint2t(vsp_velocity,vsp_depths,vsp_depths);

vsp_vrms = vint2vrms(vsp_velocity,vsp_times,vsp_times);
