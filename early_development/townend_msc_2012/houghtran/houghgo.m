cube=traces.data;%(1:500,1:300)/10;

mimg=cube;
nt=1000;
nr=1000;
plot_f=1;
 
[data_thrho,thvec,rho]=fhouhgh(mimg,nt,nr,plot_f);


% is hough transform of image to find lines and connect the points in
% original image
%------------------------
% inputs
% mimg  : binary input image
% nt : size of theta vector in hough columns
% nr : size of rho vector in hough rows
%  data_thrho= zeros(nr,nt) : is the output matrix of hough transform
% plot_f : plot flag
%-----------------------------
% outputs :
% data_thrho : matrix in hough spcae
% mimg2 : image after points connected to lines that was above threshold
% npthresh

imge=cube;
%map=data_thrho;
%theta=thvec;
minLen=10;
maxLineNum=100;
maxPeakNum=100;
resol=100;
maxGap=10;

bw=cube;
mapSz=[1000,1000];

[map,theta,rho]=hough1(bw,mapSz);
%HOUGH1 Do the Hough transform
%   [MAP,THETA,RHO] = HOUGH1(BW,MAPSZ)	BW is the input binary image, MAPSZ
%   is the size of the output MAP, means that the number of bins of rho and
%   theta is MAPSZ(1) and MAPSZ(2). Output MAP(I,J) is the number of points
%   which satisfy RHO(I) = x*sin(THETA(J))+y*cos(THETA(J)), the origin is the
%   bottom-left corner, THETA is in [-90,90).
%		If MAPSZ = [], it will be calculated automatically.


[lines]=findlines(imge,map,theta,rho,minLen,maxLineNum,maxPeakNum,resol,maxGap);
%FINDLINE Find line features in the image using Hough transform
%   LINES = FINDLINES(IMGE,MAP,THETA,RHO,MINLEN,MAXLINENUM,MAXPEAKNUM,RESOL,MAXGAP)	
%	IMGE is the binary edge image; 
%	MAP,THETA,RHO is the output of HOUGH;
%	MINLEN is the minimum length of lines that will be found;
%	No more MAXLINENUM lines will be returned;
%	MAXPEAKNUM is the maximum peaks	that will be detected on the map;
%	MAXGAP is the maximum gap length that will be merged;
%	RESOL is the resolution of the detected lines, the smaller, the lines
%	must be further away from each other. A suggest value is 100.
%	LINES is a struct array with members: theta,rho,point1([row,col]),point2. 


%plot lines extracted
figure('Color',[1 1 1])
for i=1:size(lines,2)
    plot([lines(i).point1(2),lines(i).point2(2)],...
        [lines(i).point1(1),lines(i).point2(1)])
    hold on
end
set(gca,'YDir','Reverse')


