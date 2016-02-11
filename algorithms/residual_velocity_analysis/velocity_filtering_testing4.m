
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
load('/segy/URY/2014_BG_water_column_imaging/matlab/gridfit_rms_picks2_800x40.mat');

% try calculating interval sample-by-sample straight from gridded rms


% figure; scatter(alllocs,-alltimes,20,allvels,'filled'); caxis([1.48 1.54]);
% figure; scatter(alllocs,-alltimes,20,allsmvels,'filled'); caxis([1.48 1.54]);

% dix interval velocities

vint_dix=zeros(500,dec_traces);

for trace = 1:dec_traces;
    
    ttimes = (10:10:5000)';
    vrms = gridtv(:,trace);
    
    tts = [0;ttimes(1:end-1)];
    vrms_s = [vrms(1);vrms(1:end-1)];
    
    vint_dix(:,trace) = sqrt(((vrms.^2.*ttimes) - (vrms_s.^2.*tts)) ./ (ttimes-tts));
    
end

% figure; scatter(alllocs,-alltimes,20,vint_dix,'filled'); caxis([1.48 1.54]);

figure;imagesc(vint_dix);caxis([1.48 1.54]);

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

