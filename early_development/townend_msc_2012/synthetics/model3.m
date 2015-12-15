function [data,cube,Lcube_true]=model3(xmax,ymax,zmax,nfaults,amax,bmax)

% set up output cube
xmin=1;
ymin=1;
zmin=1;

cube=zeros(zmax-zmin+1,xmax-xmin+1,ymax-ymin+1);
Lcube_true=cube;

for f=1:nfaults
    f
    a=3+rand*abs(amax-3);
    b=3+rand*abs(bmax-3);
    xc=rand*xmax;
    yc=rand*ymax;
    zc=rand*zmax;
    theta=45;%rand*pi;%deg2rad(70)+rand*pi/4.5;
    phi=rand*pi;
    PHI=70;%deg2rad(70)+rand*pi/4.5;
    % create ellipse at origin
    coord=[];
    for x=-a:0.5:a
        for y=-b:0.5:b
            ell=(x/a)^2+(y/b)^2;
            if ell<=1
                coord=[coord;[0,x,y]];
            end
            
        end
    end
    
    % rotate ellipse
    R=[cos(theta)*cos(PHI),-cos(phi)*sin(PHI)+sin(phi)*sin(theta)*cos(PHI),...
        sin(phi)*sin(PHI)+cos(phi)*sin(theta)*cos(PHI);...
        
        cos(theta)*sin(PHI),cos(phi)*cos(PHI)+sin(phi)*sin(theta)*sin(PHI),...
        -sin(phi)*cos(PHI)+cos(phi)*sin(theta)*sin(PHI);...
        
        -sin(theta),sin(phi)*cos(theta),cos(phi)*cos(theta)];
    
    coord=(R*coord')';
    
    % move ellipse to correct position and discretise
    coord=coord+repmat([zc,xc,yc],size(coord,1),1);
    coord=round(coord);
    
    %find coordinates within cube that are outside and remove
    for i=1:size(coord,1),
        if sum(coord(i,:)<[zmin,xmin,ymin])>0 ||...
                sum(coord(i,:)>[zmax,xmax,ymax])>0
            coord(i,:)=NaN;
        end
    end
    
    %get max height of the ellipse, within the cube
    len=min(max(coord(:,1)),zmax)-max(min(coord(:,1)),zmin)+1;
    
    % put coordinates into cube form, deleting entries outwith limits
    for i=1:size(coord,1),
        if isnan(coord(i,1))~=1
            if Lcube_true(coord(i,1),coord(i,2),coord(i,3))<len
                Lcube_true(coord(i,1),coord(i,2),coord(i,3))=len;
                cube(coord(i,1),coord(i,2),coord(i,3))=10;
            end
        end
    end
    
end

% unravel length cube into trace procession
data=zeros(zmax-zmin+1,(xmax-xmin+1)*(ymax-ymin+1));
for i=xmin:xmax
    ind=(ymin+(i-1)*ymax):1:(i*ymax);
    data(:,ind)=cube(:,i-xmin+1,:);
end
