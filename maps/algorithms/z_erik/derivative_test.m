xdata=0:25;
xdata1=0:24;
xdata2=0:23;
ydata=[0,4,7,9,10,9,7,4,0,-4,-7,-9,-10,-9,-7,-4,0,4,7,9,10,9,7,4,0,-4];

deriv1=diff(ydata);
deriv2=diff(ydata,2);

plot(xdata,ydata)
figure
plot(xdata1,deriv1)
figure
plot(xdata2,deriv2)
















