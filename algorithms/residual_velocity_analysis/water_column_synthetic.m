
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
