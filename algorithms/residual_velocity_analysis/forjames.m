load('/data/URY/segy/2014_BG_water_column_imaging/conditioned_datasets/cdps/4150A216/gridfit_rms_picks_200x40.mat')
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
figure;imagesc(gridtv);caxis([1.48 1.54]);

figure;scatter(alllocs,-alltimes,50,allvels,'filled'); caxis([1.48 1.54]);
 
 
gridtv_less=gridfit(alllocs,alltimes,allvels,dec_ilxl(:,2),10:10:5000,'smoothness',[50 20]);


figure; imagesc(gridtv_less);  caxis([1.48 1.54]);

gridtv_lessb=gridfit(alllocs,alltimes,allvels,dec_ilxl(:,2),10:10:5000,'smoothness',[50 20],'interp','bilinear');

figure; imagesc(gridtv_lessb);  caxis([1.48 1.54]);