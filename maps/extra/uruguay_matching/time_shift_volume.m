% Create time-shift volume to align BG and Total volumes

% import shifts from cross-correlations

shifts = importdata('/data/URY/opencps/URY_BG_Geom/tables/total_bg_time_shifts_5windows.csv');

%finl = min(shifts.data(:,1)); linl = max(shifts.data(:,1));
%fxl = min(shifts.data(:,2)); lxl = max(shifts.data(:,2));

finl = 5000; linl = 12200; fxl = 12000; lxl = 14500;


shifts_grid = zeros(5,(linl-finl)+1,(lxl-fxl)+1);
shifts_medfilt = zeros(5,(linl-finl)+1,(lxl-fxl)+1);
shifts_smth = zeros(5,(linl-finl)+1,(lxl-fxl)+1);

% filter the shifts
figure;

for window = 1:5
    
    shifts_grid(window,:,:) = griddata(shifts.data(:,2),shifts.data(:,1),shifts.data(:,(4+window)),fxl:lxl,(finl:linl)','nearest');
    shifts_medfilt(window,:,:) = medfilt2(squeeze(shifts_grid(window,:,:)),[51 51]);
    shifts_smth(window,:,:) = gaussian_2dsmth(squeeze(shifts_medfilt(window,:,:)),21,21);
    subplot(1,5,window); imagesc(squeeze(shifts_smth(window,:,:))); caxis([-20 20]);
end

wbtimes = importdata('/data/URY/opencps/URY_BG_Geom/tables/Waterbottom_Pick_Smth_4x4.csv');

wbfinl = min(wbtimes.data(:,1));
wblinl = max(wbtimes.data(:,1));
wbninl = (wblinl - wbfinl)/4+ 1;
wbfxl = min(wbtimes.data(:,2));
wblxl = max(wbtimes.data(:,2));

nlocs = size(shifts.data,1);
times = zeros(5,nlocs);
samples = zeros(5,nlocs);

% some locations don't have valid shifts after the smoothing etc
% so check for NaNs and only increment the counter if there are none
% if counter doesn't increment then bad shifts get overwritten by next
% location

% count = 1;
% 
% for loc = 1:nlocs
%     inl(count) = shifts.data(loc,1);
%     xl(count) = shifts.data(loc,2);
%     inl4 = round(inl(count)/4)*4;
%     xl4 = round(xl(count)/4)*4;
%     wb_idx = ((xl4 - wbfxl)/4)*wbninl+(inl4-wbfinl)/4+1;
%     wbt = wbtimes.data(wb_idx,3);
%     times(:,count)= [0;wbt+1000;wbt+1500;wbt+2500;wbt+3500];
%     
%     shifts_r = inl(count)-finl+1; shifts_c = xl(count)-fxl+1;
%     samples(:,count) = [shifts_smth(1,shifts_r,shifts_c);shifts_smth(1:4,shifts_r,shifts_c)];
%     
%     count=count+not(sum(isnan(samples(:,count))));
%     
% %     if(sum(isnan(samples))==0)
% %         count=count+1;
% %     end
%         
% end
% 
% count=count-1;

out_poly = importdata('/data/URY/opencps/URY_BG_Geom/tables/Timeshift_output.poly');

finlout = min(out_poly.data(:,2)); 
linlout = max(out_poly.data(:,2)); 
fxlout = min(out_poly.data(:,1)); 
lxlout = max(out_poly.data(:,1)); 

[inlgrid xlgrid] = meshgrid(finlout:linlout,fxlout:lxlout);

outlim = inpolygon(xlgrid,inlgrid,out_poly.data(:,1),out_poly.data(:,2));

inl_locs = inlgrid(outlim);
xl_locs = xlgrid(outline);

nlocs = size(inl_locs(:),1);

for loc = 1:nlocs
    inl = inl_locs(loc);
    xl = xl_locs(loc);
    inl4 = round(inl/4)*4;
    xl4 = round(xl/4)*4;
    wb_idx = ((xl4 - wbfxl)/4)*wbninl+(inl4-wbfinl)/4+1;
    wbt = wbtimes.data(wb_idx,3);
    times(:,loc)= [0;wbt+1000;wbt+1500;wbt+2500;wbt+3500];
    
    shifts_r = inl-finl+1; shifts_c = xl-fxl+1;
    samples(:,loc) = [shifts_smth(1,shifts_r,shifts_c);shifts_smth(1:4,shifts_r,shifts_c)];
                
end

 
traces = zeros(562,nlocs);

for trace=1:nlocs;
    traces(:,trace)=interp1(times(:,trace),samples(:,trace),([16:16:8992])','linear','extrap');
end

% write the output using super-easy and intuitive node_segy_write

output_vol{1,1} = 'Meta data for output files';
output_vol{2,2} = traces; % store output traces in inline/xline order
output_vol{1,2}{1,1} = [inlgrid(:),xlgrid(:)]; % store inline/xline values
output_vol{1,2}{2,1} = uint32(zeros(size(inlgrid(:)),1)); % what is this for? offsets?
output_vol{1,1} = strcat('timeshift traces ',date);
output_vol{1,3} = 'is_gather'; % 1 is yes, 0 is no
output_vol{2,3}=0;
output_vol{2,1} = 'timeshift_v3';

output_dir = '/data/URY/segy/2013_pgs_uruguay_processing/BG_Total_Merge_Working/';
% check to see if the directory exists
if exist(output_dir,'dir') == 0
    mkdir(output_dir);
end



node_segy_write(output_vol,'1',16,output_dir); % write the output
