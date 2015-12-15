%track fault identified in main
function [IDcube,fault]=track(cube,locy,locx,ID,IDcube,fault)

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
    if cube(locy,locx)==10 && IDcube(locy,locx)==0
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
while on==1
    if locy>1
        locy=locy-1;
        [IDcube locy xmin xmax on fault]=step(cube,locy,xmin,xmax,ID,IDcube,fault);
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
function [IDcube,locy,xmin_up,xmax_up,on,fault]=step(cube,locy,xmin,xmax,ID,IDcube,fault)

on=0;
xmin_up=999999;
xmax_up=0;
%check within step up range (xmin -> max)
for locx=xmin:xmax
    if cube(locy,locx)==10 && IDcube(locy,locx)==0
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

%check outwith step up range
%in the decreasing x direction
outwith=1;
locx=xmin_up+1;
while outwith==1 && on==1
    locx=locx-1;
    if locx<=0
        outwith=0;
    elseif cube(locy,locx)==10 && IDcube(locy,locx)==0
        IDcube(locy,locx)=ID;
        fault(ID).loc=[fault(ID).loc;[locx,locy]];
        xmin_up=locx-1;
    else
        outwith=0;
    end
end

%in the increasing x direction
outwith=1;
locx=xmax_up-1;
while outwith==1 && on==1
    locx=locx+1;
    if locx>size(cube,2)
        outwith=0;
    elseif cube(locy,locx)==10 && IDcube(locy,locx)==0
        IDcube(locy,locx)=ID;
        fault(ID).loc=[fault(ID).loc;[locx,locy]];
        xmax_up=locx+1;
    else
        outwith=0;
    end
end
