edge_one = importdata('/data/URY/opencps/URY_BG_Geom/tables/Timeshift_Full_Edge.poly');

edge_zero =  importdata('/data/URY/opencps/URY_BG_Geom/tables/Timeshift_Zero_Edge.poly');

edge_one.data(end+1,:)=edge_one.data(1,:);
edge_zero.data(end+1,:)=edge_zero.data(1,:);

% note full edge polygon is inline,crossline and zero edge is
% crossline,inline

count = 1;

for segment = 1:size(edge_one.data,1)-1
    xstart = edge_one.data(segment,2);
    xend = edge_one.data(segment+1,2);    
    ystart = edge_one.data(segment,1);
    yend = edge_one.data(segment+1,1);
    
    xinc = 8*((xend<xstart)*-2+1);
    num_pts = floor((xend-xstart)/xinc);

    edge_one_interp(count,:) = [xstart,ystart];
      
    x = xstart+xinc;
    
    for pt = 1:num_pts
        count=count+1;
       
        y=((x-xstart)/(xend-xstart))*(yend-ystart)+ystart;
        edge_one_interp(count,:)=[x,y];
        x = x + xinc;
        
    end
    %scatter(edge_one_interp(:,1),edge_one_interp(:,2));
end
figure; scatter(edge_one_interp(:,1),edge_one_interp(:,2));
count = 1;

for segment = 1:size(edge_zero.data,1)-1
    xstart = edge_zero.data(segment,1);
    xend = edge_zero.data(segment+1,1);    
    ystart = edge_zero.data(segment,2);
    yend = edge_zero.data(segment+1,2);
    
    ysign = (yend<ystart)*-2+1;
    xinc = 8*((xend<xstart)*-2+1);
    num_pts = floor((xend-xstart)/xinc);

    edge_zero_interp(count,:) = [xstart,ystart];
      
    x = xstart+xinc;
    
    for pt = 1:num_pts
        count=count+1;
       
        y=((x-xstart)/(xend-xstart))*(yend-ystart)+ystart;
        edge_zero_interp(count,:)=[x,y];
        x = x + xinc;
        
    end
    %scatter(edge_one_interp(:,1),edge_one_interp(:,2));
end
hold all; scatter(edge_zero_interp(:,1),edge_zero_interp(:,2));

edge_one_interp(:,3)=1;
edge_zero_interp(:,3)=0;

gridpoints = [edge_one_interp;edge_zero_interp];
xmin=floor(min(gridpoints(:,1))/8)*8;
xmax=ceil(max(gridpoints(:,1))/8)*8;
ymin=floor(min(gridpoints(:,2))/8)*8;
ymax=ceil(max(gridpoints(:,2))/8)*8;

tapergrid = griddata(gridpoints(:,1),gridpoints(:,2),gridpoints(:,3),[xmin:2:xmax],([ymin:2:ymax])');

[taperx,tapery] = meshgrid(xmin:2:xmax,ymin:2:ymax);

dlmwrite('/data/URY/opencps/URY_BG_Geom/tables/Timeshift_Taper_Grid_xlil.csv',[taperx(:),tapery(:),tapergrid(:)]);
