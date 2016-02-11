
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
% get list of lines

% linedir = '/data/URY/segy/2014_BG_water_column_imaging/conditioned_datasets/semblance/';
% filescan = 'gath_block1_semblance_io_s1c5_';
% 
% [infiles, nfiles] = directory_scan({linedir},filescan);
% 
count=0;
time1 = 2000;
time2 = 3000;
dist1 = time1*1.5;
dist2 = time2*1.5;

gridfit_hsmth = 200;
gridfit_vsmth = 20;
hsmth=801;
vsmth=41;
% 
% for file = 1:nfiles
%     if ~isempty(strfind(infiles.names{file},'.segy'))
%         count=count+1;  
%         linenames{count} = regexprep(infiles.names{file},'gath_block1_semblance_io_s1c5_','');
%         linenames{count} = regexprep(linenames{count},'.segy','');
% 
%     end
% end
% 
% nlines=size(linenames,2);

linename = '3590A159';
linedir = '/data/URY/segy/2014_BG_water_column_imaging/conditioned_datasets/swath_semblance/';
% linedir = '/data/URY/segy/2014_BG_water_column_imaging/conditioned_datasets/semblance/';

% for line=1:nlines;
    % line=45;
    disp([num2str(line),'...',linename,' .']);
    load(strcat(linedir,linename,'/semblance_picks_',linename,'.mat'));
    
    % scale from m/s to m/ms
    allvels = allvels./1000;
    
    % omit low fold locations
    fold = importdata(strcat('/data/URY/segy/2014_BG_water_column_imaging/textfiles/navseis_rad_tfd_kfilt_cdps_fold_',linename,'.csv'));
    locs = fold.data(find(fold.data(:,6)>50),1);
    fflocs=alllocs(ismember(alllocs,locs));
    ffvels = allvels(ismember(alllocs,locs));
    fftimes=alltimes(ismember(alllocs,locs));
    fcdp = round(locs(1)/10)*10;
    lcdp = round(locs(end)/10)*10;
    grididx=fftimes<=(time2/2);
    
    % import direct arrival vels
    direct = importdata(strcat('/data/URY/opencps/blocks_8_9_13/tables/direct_arrival_vel_',linename,'.txt'));
    directlocs = direct.data(:,2);
    directvels = direct.data(:,3);
    ffdirectlocs = directlocs(ismember(directlocs,locs));
    ffdirectvels = directvels(ismember(directlocs,locs));
    
    gridlocs=[fflocs(grididx);ffdirectlocs];
    gridtimes=[fftimes(grididx);repmat(50,size(ffdirectlocs,1),1)];
    gridvels=[ffvels(grididx);ffdirectvels];
    rmsgrid=gridfit(gridlocs,gridtimes,gridvels,fcdp:10:lcdp,10:10:(time2/2),'smoothness',[gridfit_hsmth gridfit_vsmth]);
    smthrms=gaussian_2dsmth(rmsgrid,hsmth,vsmth);
    vintgrid=zeros(size(smthrms));
    vintgrid=bsxfun(@vrms2vint,smthrms,(10:10:(time2/2))');
        
    wbstat1 = dist1.*(1.5^-1-(smthrms(time1/20,:)).^-1);
    wbstat2 = dist2.*(1.5^-1-(smthrms(time2/20,:)).^-1);
    
    save(strcat(linedir,linename,'/velocity_grids_',num2str(hsmth),'x',num2str(vsmth),'_',...
        num2str(time2),'ms_',linename,'.mat'),'fcdp','lcdp','rmsgrid','smthrms','vintgrid','wbstat1','wbstat2');
    
    wbstatic = [repmat(str2num(linename(6:8)),(lcdp-fcdp)/10+1,1) (fcdp:10:lcdp)' ...
        rmsgrid(time1/20,:)' rmsgrid(time2/20,:)' wbstat1' wbstat2'];
    csvwrite(strcat(linedir,linename,'/rms_',num2str(hsmth),'x',num2str(vsmth),'_',linename,'.csv'),wbstatic);
    clear locs; wbstatic; fflocs; ffvels; fftimes; fold;
%end






