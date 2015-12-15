function fault=sew_complex(fault,err_limit,var_limit,siz_limit)

%warning('off','msgid')

on=1;
match=0;
%iterate matching until no more matches can be made
while on==1
    on=0;
    ID=0;
    while ID<size(fault,2)
        ID=ID+1;
        IDcheck=ID;
        while IDcheck<size(fault,2)
            IDcheck=IDcheck+1;
            tmp_loc=[fault(ID).loc;fault(IDcheck).loc];
            if var(fault(ID).loc(:,2))>var_limit && var(fault(IDcheck).loc(:,2))>var_limit...
                    && size(fault(ID).loc,1)>siz_limit && size(fault(IDcheck).loc,1)>siz_limit
                tmp_line=polyfit(tmp_loc(:,2),tmp_loc(:,1),1);
                tmp_error=0;
                for i=1:size(tmp_loc,1)
                    tmp_offset=(tmp_line(1)*tmp_loc(i,2)+tmp_line(2)-tmp_loc(i,1))^2;
                    tmp_error=tmp_error+tmp_offset/size(tmp_loc,1);
                end
                %check requirements to match and merge if met
                if tmp_error<err_limit
                    fault(ID).loc=[fault(ID).loc;fault(IDcheck).loc];
                    fault(ID).ID=[fault(ID).ID,fault(IDcheck).ID];
                    fault(IDcheck)=[];
                    on=1;
                    match=match+1
                    ends=[0,tmp_line(2);100,(tmp_line(1)*100+tmp_line(2))];
                    plot(ends(:,2),ends(:,1));
                    hold on
                end
            end
        end
    end
end

%warning('on','msgid')