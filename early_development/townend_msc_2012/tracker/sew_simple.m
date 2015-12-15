function fault=sew_simple(fault,del_theta,del_loc,var_limit)

%warning('off','msgid')

for ID=1:size(fault,2)
    if var(fault(ID).loc(:,2))>var_limit
        %get equation of straight line fit to fault
        fault(ID).line=polyfit(fault(ID).loc(:,2),fault(ID).loc(:,1),1);
        %sort fault location data by position (useful for matching ends)
        [tmp ind]=sort(fault(ID).loc,1);
        fault(ID).loc=fault(ID).loc(flipud(ind(:,2)),:);
    else
        fault(ID).line=[inf,inf];
    end
end

on=1;
match=0;
%iterate matching until no more matches can be made
while on==1
    on=0;
    ID=0;
    while ID<size(fault,2)
        ID=ID+1;
        IDcheck=1;
        while IDcheck<=size(fault,2) && ID<=size(fault,2)
            if IDcheck~=ID
                %get difference in gradients
                diff_theta=abs(atan(fault(ID).line(1))/pi*180-...
                    atan(fault(IDcheck).line(1))/pi*180);
                %get difference between start and end points
                topbot=sqrt((fault(ID).loc(end,1)-fault(IDcheck).loc(1,1))^2+...
                    (fault(ID).loc(end,2)-fault(IDcheck).loc(1,2))^2);
                bottop=sqrt((fault(ID).loc(1,1)-fault(IDcheck).loc(end,1))^2+...
                    (fault(ID).loc(1,2)-fault(IDcheck).loc(end,2))^2);
                diff_loc=min(topbot,bottop);
                %check requirements to match and merge if met
                if diff_theta<=del_theta && diff_loc<=del_loc ...
                        && fault(IDcheck).line(1)~=Inf && fault(ID).line(1)~=Inf
                    fault(ID).loc=[fault(ID).loc;fault(IDcheck).loc];
                    fault(ID).ID=[fault(ID).ID,fault(IDcheck).ID];
                    on=1;
                    match=match+1;
                    %update line of best fit
                    fault(ID).line=polyfit(fault(ID).loc(:,2),fault(ID).loc(:,1),1);
                    %re-sort fault location data
                    [tmp ind]=sort(fault(ID).loc,1);
                    fault(ID).loc=fault(ID).loc(flipud(ind(:,2)),:);
                    %remove segment which has been merged to ID
                    fault(IDcheck)=[];
                else
                    IDcheck=IDcheck+1;
                end
            else
                IDcheck=IDcheck+1;
            end
        end
    end
end

%warning('on','msgid')