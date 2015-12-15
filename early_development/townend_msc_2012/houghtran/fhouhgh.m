% fhough=inline('x*cos(th)+y*sin(th)','x','y','th');
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
function [data_thrho,thvec,rho]=fhouhgh(mimg,nt,nr,plot_f)

data_thrho=[];mimg2=[];
%% start hough transform
% find points
ind_1 = find(mimg);np=length(ind_1);
[ii,jj]=ind2sub(size(mimg),ind_1);
% hough function
fhough=inline('x*cos(th)+y*sin(th)','x','y','th');
fdeg2rad=inline('ang*pi/180','ang');
drad2deg=inline('angrad*180/pi','angrad');
data_thrho=zeros(nr,nt);%the hough matrix initilized

% build vectors in space theta/rho of hough
thvec=fdeg2rad(linspace(-90,90,nt));
[mr,nc]=size(mimg);
rr=(sqrt(mr.^2+nc.^2));
rho=linspace(-rr,rr,nr);

%% build hough matrix using transformation
for ind=1:nt
    cur_th=thvec(ind);
    % hough
    g=fhough(jj,ii,cur_th);
    % nearest neighboor 1d
     [ind_curg]=nearest_1(rho,g);
     % cell increment
     for icell=1:length(ind_curg)
         data_thrho(ind_curg(icell),ind)=data_thrho(ind_curg(icell),ind)+1;
     end
end

if(plot_f==1)
 figure('Color',[1 1 1]),imagesc(drad2deg(thvec),rho,(data_thrho)); axis on;grid on;colorbar;
 xlabel('theta'),ylabel('rho');title('hough transform')
end


function [ind_curg]=nearest_1(rho,g)
%find index of [rho] which holds value nearest to each of [g]

ind_curg=zeros(size(g));
for i=1:size(g,1)
    difsmall=999999;
    for j=1:size(rho,2)
        dif=abs(g(i)-rho(j));
        if dif<difsmall
            ind_curg(i)=j;
            difsmall=dif;
        end
    end
end