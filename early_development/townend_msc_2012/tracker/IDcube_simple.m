%  clear fault
% % %cube=traces.data(1:500,1:500);
%  thres_len=3;

function Lcube=IDcube_simple(cube,thres_len)

%create fault ID cube
IDcube=zeros(size(cube));
ID=0;
fault.loc=[];
fault.len=[];
for i=size(cube,1):-1:1
    for j=1:size(cube,2)
        if cube(i,j) == 10 && IDcube(i,j)==0
            ID=ID+1;
            fault(ID).ID=ID;
            [IDcube,fault]=track_wfit(cube,i,j,ID,IDcube,fault);
            fault(ID).len=fault(ID).loc(1,2)-fault(ID).loc(end,2);
        end
    end
end

% remove fault segments below threshold and sort fault segments by length
ind=find([fault.len]<thres_len);
for i=1:size(ind,2)
    IDcube(find(IDcube==fault(ind(i)).ID))=0;
end
fault(ind)=[];
[tmp ind]=sort([fault.len]);
fault=fault(fliplr(ind));

%get LINEST fit to each fault, sew together faults that are similar orientation 
%and close, get fault length and height
fID=zeros(size(fault,2),1);

del_theta=5;    %bounds on sewing faults
del_loc=10;
var_limit=0.1;  %threshold variance to recoginse fault segment
err_limit=10;
siz_limit=20;
fault=sew_simple(fault,del_theta,del_loc,var_limit);
%fault=sew_complex(fault,err_limit,var_limit,siz_limit);

for ID=1:size(fault,2)
    fault(ID).len=max(fault(ID).loc(:,2))-min(fault(ID).loc(:,2))+1;
    for i=1:size(fault(ID).ID,2)
        fID(fault(ID).ID(i))=ID;
    end
end


%create fault length cube and IDcube2 which holds fault struct IDs for sewn
%faults
Lcube=zeros(size(IDcube));
%IDcube2=zeros(size(IDcube));
for i=1:size(IDcube,1)
    for j=1:size(IDcube,2)
        if IDcube(i,j)~=0 
            Lcube(i,j)=fault(fID(IDcube(i,j))).len;
            IDcube2(i,j)=fID(IDcube(i,j));
        end
    end
end

%display
% figure('Color',[1 1 1]);
% imagesc(IDcube)
% colorbar

% figure('Color',[1 1 1]);
% imagesc(-IDcube2)
% colorbar
% 
figure('Color',[1 1 1]);
imagesc(squeeze(Lcube(100,:,:)))
daspect([1 1 1])
c=colorbar;
ylabel(c,'Vertical extent (samples)');
xlabel(gca,'In-line number');
ylabel(gca,'Cross-line number');
% print('-dpng','-r300','testing2_tracker_fitting_distance.png')