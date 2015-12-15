function [datai] = it2di(data,range_in,range_out)
% Interpolates a 3D array using iterative 1D interpolation. range_in is a
% cell array containing three rows which hold the row, column and page
% vectors at which data exists. range_out is a cell array containing three
% rows which hold the row, column and page vectors onto which the data will
% be interpolated. data is a 3D array, which has length dimensions eqaul to
% the length of the vectors held in range_in.

x = range_in{1};
y = range_in{2};
z = range_in{3};

xi = range_out{1};
yi = range_out{2};
zi = range_out{3};

xi_tmp = x;
yi_tmp = y;
zi_tmp = z;

x_tmp = x;
y_tmp = y;
z_tmp = z;

f = 1;
while length(xi_tmp) <= length(x)
    xscale = 2^(floor(log2(length(xi)))-f);
    xi_tmp = downsample(xi,xscale);
    f = f+1;
end
g = 1;
while length(yi_tmp) <= length(y)
    yscale = 2^(floor(log2(length(yi)))-g);
    yi_tmp = downsample(yi,yscale);
    g = g+1;
end
h = 1;
while length(zi_tmp) <= length(z)
    zscale = 2^(floor(log2(length(zi)))-h);
    zi_tmp = downsample(zi,zscale);
    h = h+1;
end

datai0 = data;
datai1 = zeros(length(zi_tmp),length(y_tmp),length(x_tmp));
datai = zeros(length(zi_tmp),length(yi_tmp),length(xi_tmp));
counter = 1;
xskip = 0;
yskip = 0;
zskip = 0;

while xscale > 0
    for i = 1:length(y_tmp)
        for j = 1:length(x_tmp)
            if counter <= 1
                datai1(:,i,j) = interp1(z_tmp,datai0(:,i,j),zi_tmp,'nearest','extrap');
            elseif zskip == 1
                datai1(:,i,j) = datai0(:,i,j);
            else
                datai1(:,i,j) = interp1(z_tmp,datai0(:,i,j),zi_tmp,'pchip','extrap');
            end
        end
    end
    for i = 1:length(zi_tmp)
        datai1_zslice = reshape(datai1(i,:,:),length(y_tmp),length(x_tmp));
        if counter <= 1
            datai(i,:,:) = interp2(x_tmp,y_tmp,datai1_zslice,xi_tmp',yi_tmp,'nearest',mean(mean(datai1_zslice)));
        elseif ((xskip+yskip) == 2)
            datai(i,:,:) = datai1_zslice;
        else 
            datai(i,:,:) = interp2(x_tmp,y_tmp,datai1_zslice,xi_tmp',yi_tmp,'spline',mean(mean(datai1_zslice)));
        end
    end
    fprintf('Interpolation to scale %d complete in inlne direction\n',xscale);
    fprintf('Interpolation to scale %d complete in crossline direction\n',yscale);
    fprintf('Interpolation to scale %d complete in Z direction\n\n',zscale);
    if xscale+yscale+zscale == 3
        fprintf('Interpolation to scale 1 complete in all directions\n');
        break
    end
    if (((xscale+yscale) == 2) || ((xscale+zscale) == 2) || (((yscale + zscale)) == 2))
        xscale = 2;
        yscale = 2;
        zscale = 2;
    end
    if xscale~=1
        xscale = xscale/2;
    else
        xskip = 1;
    end
    if yscale~=1
        yscale = yscale/2;
    else
        yskip = 1;
    end
    if zscale~=1
        zscale = zscale/2;
    else
        zskip = 1;
    end
    
    x_tmp = xi_tmp;
    xi_tmp = downsample(xi,xscale);
    y_tmp = yi_tmp;
    yi_tmp = downsample(yi,yscale);
    z_tmp = zi_tmp;
    zi_tmp = downsample(zi,zscale);
    datai0 = datai;
    datai1 = zeros(length(zi_tmp),length(y_tmp),length(x_tmp));
    datai = zeros(length(zi_tmp),length(yi_tmp),length(xi_tmp));
    counter = counter+1;
end
end