% make smooth waterbottom map for keying matching filters

% should have kept everything on 4x4 grid but left 1x1 in griddata without thinking

% import picks
wbt = importdata('/data/URY/opencps/URY_BG_Geom/tables/Export_Waterbottom_Pick_4x4.csv');

% check inl/xl limits
wbfxl = min(wbt.data(:,2));
wblxl = max(wbt.data(:,2));
wbfinl = min(wbt.data(:,1));
wblinl = max(wbt.data(:,1));

% import survey outline
outline=importdata('/data/URY/opencps/URY_BG_Geom/tables/BG_Outline.poly');

% limit input data to survey outline
wbt_lim = inpolygon(wbt.data(:,2),wbt.data(:,1),outline.data(:,1),outline.data(:,2));
wbt_poly = [wbt.data(wbt_lim,2),wbt.data(wbt_lim,1),wbt.data(wbt_lim,5)];

% grid input data to form matrix for filtering
wbtgrid2 = griddata(wbt_poly(wbt_poly(:,3)<5000,1),wbt_poly(wbt_poly(:,3)<5000,2),wbt_poly(wbt_poly(:,3)<5000,3),wbfxl:wblxl,(wbfinl:wblinl)','nearest');

% median filter to remove bad picks
wbtgrid2_med101 = medfilt2(wbtgrid2,[101 101]);

% smoothing to remove any sharp changes
wbtgrid2_med101_smth101 = gaussian_2dsmth(wbtgrid2_med101,101,101);

% limit smoothed grid to outline - this requires making x,y,z vectors from
% the grid first
[r c] = find(wbtgrid);
outline_r = outline.data(:,2)-wbfinl;
outline_c = outline.data(:,1)-wbfxl;
gridlim=inpolygon(r,c,outline_r,outline_c);
gridlim = reshape(gridlim,9401,12801);
wbtgrid_filt_lim = wbtgrid2_med101_smth101.*gridlim;

% finally extrapolate to get values for all inlines/crosslines
% need x,y,z vectors of non-zero values
[wbtgrid_filt_lim_y wbtgrid_filt_lim_x]  = find(wbtgrid_filt_lim);
wbtgrid_filt_lim_z = wbtgrid_filt_lim(wbtgrid_filt_lim>0);
wbtgrid_final = griddata(wbtgrid_filt_lim_x,wbtgrid_filt_lim_y,wbtgrid_filt_lim_z,1:12801,(1:9401)','nearest');

% make matrix to export to text file
output_grid = [r+wbfinl c+wbfxl wbtgrid_final(:)];

% too slow to output every value so limit to every 4th
outgridtmp = [output_grid(mod(output_grid(:,1),4)==0,1) output_grid(mod(output_grid(:,1),4)==0,2) output_grid(mod(output_grid(:,1),4)==0,3)];
output_grid4 = [outgridtmp(mod(outgridtmp(:,2),4)==0,1) outgridtmp(mod(outgridtmp(:,2),4)==0,2)  outgridtmp(mod(outgridtmp(:,2),4)==0,3) ];

dlmwrite('/data/URY/opencps/URY_BG_Geom/tables/Waterbottom_Pick_Smth_4x4.csv',output_grid4);
