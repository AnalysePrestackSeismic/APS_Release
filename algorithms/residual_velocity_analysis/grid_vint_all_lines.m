
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

linedir = '/data/URY/segy/2014_BG_water_column_imaging/conditioned_datasets/cdps/';
filescan = 'gath_block1_conditioned_cdps_s1c5_';

[infiles, nfiles] = directory_scan({linedir},filescan);

count=0;

for file = 1:nfiles
    if ~isempty(strfind(infiles.names{file},'.segy'))
        count=count+1;  
        linenames{count} = regexprep(infiles.names{file},'gath_block1_conditioned_cdps_s1c5_','');
        linenames{count} = regexprep(linenames{count},'.segy','');

    end
end

nlines=size(linenames,2);

for line=1:nlines;
    disp([num2str(line),'...',linenames{line},' .']);
    matfile = strcat(linedir,linenames{line},'/gridfit_rms_picks_800x40_',linenames{line},'.mat');
    load(matfile,'alllocs','alltimes','gridtv','grid_xl','datadir','stack_all');
   
    % interval vels at pick locations =========================================
    
    % decimate original pick times onto new grid
           
    alllocs10=alllocs(mod(alllocs,10)==0);
    alltimes10=alltimes(mod(alllocs,10)==0);
    numpicks10=size(alllocs10,2);
    allsmvels10=zeros(1,numpicks10);
    
    % extract smoothed vels at picked locations
    alllocs_grid = ((alllocs10-alllocs10(1))./10)+1;
    
    for ii=1:numpicks10
        allsmvels10(ii)=gridtv(ceil(alltimes10(ii)/10),alllocs_grid(ii));
        
    end
    
    % interval velocities
    
    allvint10=zeros(size(allsmvels10));
    
    for trace = unique(alllocs10);
        picks=alllocs10==trace;
        ttimes=alltimes10(picks);
        vrms=allsmvels10(picks);
        tts = [0 ttimes(1:end-1)];
        vrms_s = [vrms(1) vrms(1:end-1)];
        allvint10(picks) = sqrt(((vrms.^2.*ttimes) - (vrms_s.^2.*tts)) ./ (ttimes-tts));
        
    end
    
    % copy each pick to immediately after previous pick
    % then interp1d to fill between
    % and reshape to make into a grid
    
    
    alltimes10_2 = zeros(size(alltimes10,2)*2-1,1);
    allvels10_2 = zeros(size(allvint10,2)*2-1,1);
    
    alltimes10_2(1:2:end)=round(alltimes10+(alllocs_grid-1).*5100);
    alltimes10_2(2:2:end-1)=round(alltimes10(1:end-1)+1+(alllocs_grid(1:end-1)-1).*5100);
    
    allvels10_2(1:2:end)=allvint10;
    allvels10_2(2:2:end-1)=allvint10(2:end);
    
    dec_traces = size(dec_xl,1);
    
    grid_vint = reshape(interp1(alltimes10_2,allvels10_2,1:5100*dec_traces),5100,dec_traces);
    
    filename = strcat(linedir,linenames{line},'/gridfit_rms_picks_800x40_vint_',linenames{line},'.mat');
    
    save(filename);
    
    % ========================================================================================
    % plot vint and stack

%     scrsz = get(0,'ScreenSize'); % get size of screen to use for controlling size of plot
%     
%     vint_plot = figure('Name',['Interval velocities line ', linenames{line}],'OuterPosition',[1 scrsz(4) scrsz(3)/2 scrsz(4)]); 
%     imagesc(grid_vint); caxis([1.48 1.54]); colorbar; hold all; 
%     scatter(alllocs_grid,alltimes10,10,'black','filled');
%     plotfile=strcat(datadir,'/interval_velocity_',linenames{line},'.png');
%     saveas(vint_plot,plotfile);
%     close(vint_plot);
%     
%     stack_plot = figure('Name',['Constant velocity stack line ', linenames{line}],'OuterPosition',[1 scrsz(4) scrsz(3)/2 scrsz(4)]);
%     imagesc(stack_all); colormap(gray); caxis([-0.1 0.1])
%     plotfile=strcat(datadir,'/const_vel_stack_',linenames{line},'.png');
%     saveas(stack_plot,plotfile);
%     close(stack_plot);
     
    
end


% gridded dix interval velocities for shallow =============================
% this is necessary to merge in the direct arrival velocities


% vint_shallow=zeros(50,dec_traces);
% 
% for trace = 1:dec_traces;
%     
%     ttimes = (10:10:500)';
%     vrms = gridtv(1:50,trace);
%     
%     tts = [0;ttimes(1:end-1)];
%     vrms_s = [vrms(1);vrms(1:end-1)];
%     
%     vint_shallow(:,trace) = sqrt(((vrms.^2.*ttimes) - (vrms_s.^2.*tts)) ./ (ttimes-tts));
%     
%     vint_merge(:,trace) = interp1([10:10:200 250:5000],[vint_shallow(1:20,trace)' grid_vint(250:5000,trace)'],1:5000,'linear','extrap');
% end






