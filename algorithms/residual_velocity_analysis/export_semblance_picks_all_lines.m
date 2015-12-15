
% get list of lines

linelist = importdata('/data/URY/segy/2014_BG_water_column_imaging/conditioned_datasets/semblance/lines_to_process.txt');

nlines = size(linelist,1);
basedir = '/data/URY/segy/2014_BG_water_column_imaging/conditioned_datasets/semblance/';

for line=1:nlines
    infile=strcat(basedir,linelist{line},'/velocity_grids_2000ms_',linenames{line},'.mat');
    load(infile);
    



linedir = '/data/URY/segy/2014_BG_water_column_imaging/conditioned_datasets/semblance/';
filescan = 'gath_block1_semblance_io_s1c5_';

[infiles, nfiles] = directory_scan({linedir},filescan);

count=0;

for file = 1:nfiles
    if ~isempty(strfind(infiles.names{file},'.segy'))
        count=count+1;  
        linenames{count} = regexprep(infiles.names{file},'gath_block1_semblance_io_s1c5_','');
        linenames{count} = regexprep(linenames{count},'.segy','');

    end
end

nlines=size(linenames,2);

for line=1:nlines;
% line=45;
    disp([num2str(line),'...',linenames{line},' .']);
    load(strcat(linedir,linenames{line},'/semblance_picks_',linenames{line},'.mat'));
    
    % omit low fold locations
    fold = importdata(strcat('/data/URY/segy/2014_BG_water_column_imaging/textfiles/navseis_rad_tfd_kfilt_cdps_fold_',linenames{line},'.csv'));
    locs = fold.data(find(fold.data(:,6)>50),1);
    fflocs=alllocs(ismember(alllocs,locs)); 
    ffvels = allvels(ismember(alllocs,locs)); 
    fftimes=alltimes(ismember(alllocs,locs));
    fcdp = round(locs(1)/10)*10;
    lcdp = round(locs(end)/10)*10;
    grididx=fftimes<=1200;
    gridlocs=fflocs(grididx);
    gridtimes=fftimes(grididx);
    gridvels=ffvels(grididx);
    rmsgrid=gridfit(gridlocs,gridtimes,gridvels,fcdp:10:lcdp,10:10:1000,'smoothness',[1600 200]);
    vintgrid=zeros(size(rmsgrid));
    vintgrid=bsxfun(@vrms2vint,rmsgrid,(10:10:1000)');
    save(strcat(linedir,linenames{line},'/velocity_grids_2000ms_',linenames{line},'.mat'),'fcdp','lcdp','rmsgrid','vintgrid');
    wbstatic = [repmat(str2num(linenames{line}(6:8)),(lcdp-fcdp)/10+1,1) (fcdp:10:lcdp)' rmsgrid(100,:)' ((1-1500./rmsgrid(100,:)).*2000)'];
    csvwrite(strcat(linedir,linenames{line},'/rms_2000ms_',linenames{line},'.csv'),wbstatic);
    clear locs; wbstatic; fflocs; ffvels; fftimes; fold;
    
end






