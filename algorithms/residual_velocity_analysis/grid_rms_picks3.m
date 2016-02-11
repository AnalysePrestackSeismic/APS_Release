
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

% add vels from analysis of direct arrivals

direct_velfile='/URY/opencps/blocks_8_9_13/tables/direct_arrival_vel_3790A173.txt';
directvels=importdata(direct_velfile,' ',1);
directvels_interp=interp1(directvels.data(:,2),directvels.data(:,3),ilxl(:,2));

total_picks=0;
total_traces=0;

for iblock = 2:23
    matfile = strcat('/segy/URY/2014_BG_water_column_imaging/matlab/rms_picks_test_block',int2str(iblock),'.mat');
    load(matfile,'ntraces','tvpairs','stack_traces');
    total_traces = total_traces + ntraces;
    ilxl(total_traces-ntraces+1:total_traces,:)=double(stack_traces{1,2}{1,1});
    
    for trace = 1:ntraces
        % put in picks at 1 and 100 from direct arrival vels
        % copy last pick to end for interval vel conversion
        current_trace = trace+total_traces-ntraces;
        numpicks(current_trace)=size(tvpairs{trace},1)+3; % add 3 to numpicks for direct vels and last pick padding
        pick_idx(current_trace)=total_picks+1;
        total_picks=total_picks+numpicks(current_trace);
        
        alltimes(pick_idx(current_trace):pick_idx(current_trace)+1)=[1,100];
        allvels(pick_idx(current_trace):pick_idx(current_trace)+1)=directvels_interp(current_trace);
              
        alltimes(total_picks)=5000;
        allvels(total_picks)=tvpairs{trace}(end,2);
        
        alltimes(pick_idx(current_trace)+2:total_picks-1)=tvpairs{trace}(:,1);
        allvels(pick_idx(current_trace)+2:total_picks-1)=tvpairs{trace}(:,2);
        alllocs(pick_idx(current_trace):total_picks)=ilxl(current_trace,2);
    end
    
end


% figure;scatter(alllocs,-alltimes,50,allvels,'filled'); caxis([1.48 1.54]);

% interpolate/smooth onto a grid

dec_traces = floor(total_traces/10);

%dec_trace_list = 10:10:total_traces-mod(total_traces,10);
dec_ilxl = double(ilxl(10:10:end,:));

hsmth = 800 ; vsmth = 40;

gridtv=gridfit(alllocs,alltimes,allvels,dec_ilxl(:,2),10:10:5000,'smoothness',[hsmth vsmth]);

filename=strcat('/segy/URY/2014_BG_water_column_imaging/matlab/gridfit_rms_picks2_',int2str(hsmth),'x',int2str(vsmth));

save(filename);


        