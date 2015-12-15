load('/data/URY/segy/2014_BG_water_column_imaging/conditioned_datasets/cdps/4150A216/gridfit_rms_picks_200x40.mat')

figure;imagesc(gridtv);caxis([1.48 1.54]);

figure;scatter(alllocs,-alltimes,50,allvels,'filled'); caxis([1.48 1.54]);
 
 
gridtv_less=gridfit(alllocs,alltimes,allvels,dec_ilxl(:,2),10:10:5000,'smoothness',[50 20]);


figure; imagesc(gridtv_less);  caxis([1.48 1.54]);

gridtv_lessb=gridfit(alllocs,alltimes,allvels,dec_ilxl(:,2),10:10:5000,'smoothness',[50 20],'interp','bilinear');

figure; imagesc(gridtv_lessb);  caxis([1.48 1.54]);