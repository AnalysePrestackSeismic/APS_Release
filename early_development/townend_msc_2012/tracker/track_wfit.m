%track fault identified in main
function [IDcube,fault]=track_wfit(cube,locy,locx,ID,IDcube,fault)

fit_thres=20;%20; % for full
del_theta=2;%2; % for full

%deal with first row and get 2nd row x search range
on=1;
if locx==1
    xmin=locx;
else
    xmin=locx-1;
end
xmax=locx+1;
while on==1
    on=0;
    if cube(locy,locx)~=0 && IDcube(locy,locx)==0
        IDcube(locy,locx)=ID;
        fault(ID).loc=[fault(ID).loc;[locx,locy]];
        on=1;
        xmax=locx+1;
    end
    if locx<size(cube,2)
        locx=locx+1;
    else
        on=0;
        xmax=locx;
    end
end

%step up from first row with search width of 1 trace
%return to main when no fault detected within trace window up
on=1;
go=0;
fit_curr=[];
while on==1
    if locy>1
        locy=locy-1;
        [IDcube locy xmin xmax on fault fit_curr]...
            =step(cube,locy,xmin,xmax,ID,IDcube,fault,fit_thres,fit_curr,del_theta,go);
        if size(fault(ID).loc,1)>fit_thres
            fit_curr=polyfit(fault(ID).loc(end-fit_thres:end,2),...
                fault(ID).loc(end-fit_thres:end,1),1);
            go=1;
        end
    else
        on=0;
    end
    if xmax>size(cube,2)
        xmax=size(cube,2);
    end
    if xmin<1
        xmin=1;
    end
end


%called from function track, takes one step up the fault
function [IDcube,locy,xmin_up,xmax_up,on,fault,fit_curr]...
    =step(cube,locy,xmin,xmax,ID,IDcube,fault,fit_thres,fit_curr,del_theta,go)

on=0;
xmin_up=999999;
xmax_up=0;
diff_theta=del_theta;
%check within step up range (xmin -> max)
for locx=xmin:xmax
    if cube(locy,locx)~=0 && IDcube(locy,locx)==0        
        if go==1
            fit_tmp=polyfit([fault(ID).loc(end-fit_thres:end,2);locy],...
                [fault(ID).loc(end-fit_thres:end,1);locx],1);
            diff_theta=abs(atan(fit_curr(1))/pi*180-...
                atan(fit_tmp(1))/pi*180);
        end
        if diff_theta<=del_theta
            on=1;
            IDcube(locy,locx)=ID;
            fault(ID).loc=[fault(ID).loc;[locx,locy]];
            if locx<(xmin_up)
                xmin_up=locx-1;
            end
            if locx>(xmax_up-1)
                xmax_up=locx+1;
            end
        end
    end
end

%check outwith step up range
%in the decreasing x direction
outwith=1;
locx=xmin_up+1;
diff_theta=del_theta;
while outwith==1 && on==1
    locx=locx-1;
    if locx<=0
        outwith=0;
    elseif cube(locy,locx)~=0 && IDcube(locy,locx)==0
        if go==1
            fit_tmp=polyfit([fault(ID).loc(end-fit_thres:end,2);locy],...
                [fault(ID).loc(end-fit_thres:end,1);locx],1);
            diff_theta=abs(atan(fit_curr(1))/pi*180-...
                atan(fit_tmp(1))/pi*180);
        end
        if diff_theta<=del_theta
            IDcube(locy,locx)=ID;
            fault(ID).loc=[fault(ID).loc;[locx,locy]];
            xmin_up=locx-1;
        end
    else
        outwith=0;
    end
end

%in the increasing x direction
outwith=1;
locx=xmax_up-1;
diff_theta=del_theta;
while outwith==1 && on==1
    locx=locx+1;
    if locx>size(cube,2)
        outwith=0;
    elseif cube(locy,locx)~=0 && IDcube(locy,locx)==0
        if go==1
            fit_tmp=polyfit([fault(ID).loc(end-fit_thres:end,2);locy],...
                [fault(ID).loc(end-fit_thres:end,1);locx],1);
            diff_theta=abs(atan(fit_curr(1))/pi*180-...
                atan(fit_tmp(1))/pi*180);
        end
        if diff_theta<=del_theta
            IDcube(locy,locx)=ID;
            fault(ID).loc=[fault(ID).loc;[locx,locy]];
            xmax_up=locx+1;
        end
    else
        outwith=0;
    end
end
