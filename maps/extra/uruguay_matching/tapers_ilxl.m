bg_only_poly = importdata('/data/URY/opencps/URY_BG_Geom/tables/BG_Only_Poly_ilxl.csv');

total_only_poly =  importdata('/data/URY/opencps/URY_BG_Geom/tables/Total_Only_Poly_ilxl.csv');

bg_only_poly.data(end+1,:)=bg_only_poly.data(1,:);
total_only_poly.data(end+1,:)=total_only_poly.data(1,:);

count = 1;

for segment = 1:size(bg_only_poly.data,1)-1
    xstart = bg_only_poly.data(segment,4);
    xend = bg_only_poly.data(segment+1,4);    
    ystart = bg_only_poly.data(segment,3);
    yend = bg_only_poly.data(segment+1,3);
    
    xinc = 8*((xend<xstart)*-2+1);
    num_pts = floor((xend-xstart)/xinc);

    bg_interp_poly(count,:) = [bg_only_poly.data(segment,4),bg_only_poly.data(segment,3)];
      
    x = xstart+xinc;
    
    for pt = 1:num_pts
        count=count+1;
       
        y=((x-xstart)/(xend-xstart))*(yend-ystart)+ystart;
        bg_interp_poly(count,:)=[x,y];
        x = x + xinc;
        
    end
    %scatter(bg_interp_poly(:,1),bg_interp_poly(:,2));
end
figure; scatter(bg_interp_poly(:,1),bg_interp_poly(:,2));
count = 1;

for segment = 1:size(total_only_poly.data,1)-1
    xstart = total_only_poly.data(segment,4);
    xend = total_only_poly.data(segment+1,4);    
    ystart = total_only_poly.data(segment,3);
    yend = total_only_poly.data(segment+1,3);
    
    ysign = (yend<ystart)*-2+1;
    xinc = 8*((xend<xstart)*-2+1);
    num_pts = floor((xend-xstart)/xinc);

    total_interp_poly(count,:) = [total_only_poly.data(segment,4),total_only_poly.data(segment,3)];
      
    x = xstart+xinc;
    
    for pt = 1:num_pts
        count=count+1;
       
        y=((x-xstart)/(xend-xstart))*(yend-ystart)+ystart;
        total_interp_poly(count,:)=[x,y];
        x = x + xinc;
        
    end
    %scatter(bg_interp_poly(:,1),bg_interp_poly(:,2));
end
figure; scatter(total_interp_poly(:,1),total_interp_poly(:,2));

bg_interp_poly(:,3)=1;
total_interp_poly(:,3)=0;

gridpoints = [bg_interp_poly;total_interp_poly];
xmin=floor(min(gridpoints(:,1))/8)*8;
xmax=ceil(max(gridpoints(:,1))/8)*8;
ymin=floor(min(gridpoints(:,2))/8)*8;
ymax=ceil(max(gridpoints(:,2))/8)*8;

tapergrid = griddata(gridpoints(:,1),gridpoints(:,2),gridpoints(:,3),[xmin:xmax],([ymin:ymax])');

[taperx,tapery] = meshgrid(xmin:xmax,ymin:ymax);

dlmwrite('/data/URY/opencps/URY_BG_Geom/tables/Taper_Grid_xlil.csv',[taperx(:),tapery(:),tapergrid(:)]);
