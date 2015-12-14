
% residual_velocity_analysis('/URY/segy/2014_BG_water_column_imaging/matlabjob_meta/job_meta_23Jul2014.mat','10','velocity_path');

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

line=20;

% for line=1:nlines;
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
        
%     end
    
    
    % figure;scatter(alllocs,-alltimes,50,allvels,'filled'); caxis([1.48 1.54]);
    
    % interpolate/smooth onto a grid
    
    % gridtv=gridfit(alllocs,alltimes,allvels,10:10:ceil(total_traces/10),10:10:5000,'smoothness',100);
    
    % dec_ilxl = double(ilxl(10:10:end,:));
    hsmth = 20000 ; vsmth = 100;
    gridtv=gridfit(alllocs,alltimes,allvels,ilxl(:,2),10:10:5000,'smoothness',[hsmth vsmth]);
    
    filename=strcat(datadir,'gridfit_rms_picks_',int2str(hsmth),'x',int2str(vsmth),'_',linenames{line},'.mat');
    
    save(filename);

    
end
        




        