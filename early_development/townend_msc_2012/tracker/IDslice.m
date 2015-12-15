% clear fault
% slice=traces.data(1:500,1:500);
% thres_len=3;

function Lslice=IDslice(slice)

%create fault ID slice
IDslice=zeros(size(slice));
ID=0;
fault.loc=[];
fault.len=[];
for i=size(slice,1):-1:1
    for j=1:size(slice,2)
        if slice(i,j) ~= 0 && IDslice(i,j)==0
            ID=ID+1;
            fault(ID).ID=ID;
            [IDslice,fault]=track_wfit(slice,i,j,ID,IDslice,fault);
            fault(ID).len=fault(ID).loc(1,2)-fault(ID).loc(end,2);
        end
    end
end

%sort fault segments by length
[tmp ind]=sort([fault.len]);
fault=fault(fliplr(ind));

%get LINEST fit to each fault, sew together faults that are similar orientation 
%and close, get fault length and height
fID=zeros(size(fault,2),1);

del_theta=5;    %bounds on sewing faults
del_loc=10;
var_limit=0.1;  %threshold variance to recoginse fault segment
%err_limit=10;
%siz_limit=20;
fault=sew_simple(fault,del_theta,del_loc,var_limit);
%fault=sew_complex(fault,err_limit,var_limit,siz_limit);

for ID=1:size(fault,2)
    fault(ID).val=zeros(size(fault(ID).loc,1),1);
    for i=1:size(fault(ID).loc,1)
        fault(ID).val(i)=slice(fault(ID).loc(i,2),fault(ID).loc(i,1));
    end
    fault(ID).med=max(fault(ID).val);
    for i=1:size(fault(ID).ID,2)
        fID(fault(ID).ID(i))=ID;
    end
end


%create fault length slice and IDslice2 which holds fault struct IDs for sewn
%faults
Lslice=zeros(size(IDslice));
%IDslice2=zeros(size(IDslice));
for i=1:size(IDslice,1)
    for j=1:size(IDslice,2)
        if IDslice(i,j)~=0 
            Lslice(i,j)=fault(fID(IDslice(i,j))).med;
            IDslice2(i,j)=fID(IDslice(i,j));
        end
    end
end

%display
% figure('Color',[1 1 1]);
% imagesc(IDslice)
% colorbar

% figure('Color',[1 1 1]);
% imagesc(-IDslice2)
% colorbar
% 
% figure('Color',[1 1 1]);
% imagesc((Lslice))
% colorbar