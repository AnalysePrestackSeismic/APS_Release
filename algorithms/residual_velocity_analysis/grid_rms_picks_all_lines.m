
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

% 'velocity_path' parameter is not used and is only there for future
% implementation

% grid using "gridfit" algorithm which has inherent smoothing
% then only output on 10x10 grid 

% get list of lines
datadir = '/data/URY/segy/2014_BG_water_column_imaging/conditioned_datasets/cdps/';
filescan = 'gath_block1_conditioned_cdps_s1c5_';

[infiles, nfiles] = directory_scan({datadir},filescan);

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
    clear blocknum ilxl stack_all stack_mask_all numpicks pick_idx alltimes allvels alllocs;
    disp(['Working on ', linenames{line}]);
    datadir = strcat('/data/URY/segy/2014_BG_water_column_imaging/conditioned_datasets/cdps/',linenames{line},'/');
    [pickfiles, npickfiles] = directory_scan({datadir},'rms_picks_block'); % get list of mat files
    for iblock=1:npickfiles % loop round a get block number of from each filename
        blocknum(iblock)=str2num(regexprep(regexprep(pickfiles.names{iblock},'rms_picks_block',''),'.mat',''));
    end
    blocknum = sort(blocknum);
    
    % loop round doing the rms gridding for each block
    
    total_picks=0;
    total_traces=0;
    
    for iblock = blocknum
        matfile = strcat(datadir,'rms_picks_block',int2str(iblock),'.mat');
        load(matfile,'ntraces','tvpairs','stack_traces','stack','stack_mask');
        total_traces = total_traces + ntraces;
        ilxl(total_traces-ntraces+1:total_traces,:)=double(stack_traces{1,2}{1,1});
        stack_all(:,total_traces-ntraces+1:total_traces) = stack;
        stack_mask_all(:,total_traces-ntraces+1:total_traces) = stack_mask;
        
        for trace = 1:ntraces
            % copy the first pick to the first time to avoid NaNs in gridding
            % copy last pick to end for interval vel conversion
            current_trace = trace+total_traces-ntraces;
            numpicks(current_trace)=size(tvpairs{trace},1)+2; % add 2 to numpicks for first/last padding
            pick_idx(current_trace)=total_picks+1;
            total_picks=total_picks+numpicks(current_trace);
            
            alltimes(pick_idx(current_trace))=1;
            allvels(pick_idx(current_trace))=tvpairs{trace}(1,2);
            
            alltimes(total_picks)=5000;
            allvels(total_picks)=tvpairs{trace}(end,2);
            
            alltimes(pick_idx(current_trace)+1:total_picks-1)=tvpairs{trace}(:,1);
            allvels(pick_idx(current_trace)+1:total_picks-1)=tvpairs{trace}(:,2);
            alllocs(pick_idx(current_trace):total_picks)=ilxl(current_trace,2);
            
        end
        
    end
    
    
    % figure;scatter(alllocs,-alltimes,50,allvels,'filled'); caxis([1.48 1.54]);
    
    % interpolate/smooth onto a grid
    
    % gridtv=gridfit(alllocs,alltimes,allvels,10:10:ceil(total_traces/10),10:10:5000,'smoothness',100);
    
    % dec_ilxl = double(ilxl(10:10:end,:));
    
    % dec_xl = ilxl(mod(ilxl(:,2),10)==0,2); % decimate to every 10th cdp
    
    grid_xl = (round(ilxl(1,2)/10)*10):10:(round(ilxl(end,2)/10)*10); % define output grid for gridfit as every 10th cdp
    
    hsmth = 800 ; vsmth = 40;
    gridtv=gridfit(alllocs,alltimes,allvels,grid_xl,10:10:5000,'smoothness',[hsmth vsmth]);
    
    filename=strcat(datadir,'gridfit_rms_picks_',int2str(hsmth),'x',int2str(vsmth),'_',linenames{line},'.mat');
    
    save(filename);

    
end
        




        