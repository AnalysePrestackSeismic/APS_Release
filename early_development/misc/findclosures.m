function [closuregrd] = findclosures(slice,contint)

[nxl nil] = size(slice);
closuregrd = zeros(size(slice));

C = contourc(slice,contint);

loopcount = 1;
ifcount = 1;
jump = 1;
nclosures = 0;

while jump < length(C)
    %Z(loopcount,1) = C(1,jump);
    if C(1,jump+1) == C(1,jump+C(2,jump)) && C(2,jump+1) == C(2,jump+C(2,jump))
        nclosures = nclosures+1;
        closures{ifcount,1} = C(:,jump:jump+C(2,jump));
        ifcount = ifcount+1;
    end
    jump = jump+C(2,jump)+1;
    loopcount = loopcount+1;
end

for i = 1:nclosures
    xmax = ceil(max(closures{i,1}(1,2:end)));
    ymax = ceil(max(closures{i,1}(2,2:end)));
    xmin = floor(min(closures{i,1}(1,2:end)));
    ymin = floor(min(closures{i,1}(2,2:end)));
    
    sqlength = max(xmax-xmin,ymax-ymin);
    
    xsample = (xmin:1:xmin+sqlength)';
    ysample = ymin:1:ymin+sqlength;
    
    xsamplesq = repmat(xsample,1,sqlength+1);
    ysamplesq = repmat(ysample,sqlength+1,1);
    
    IN = inpolygon(xsamplesq,ysamplesq,closures{i}(1,2:end)',closures{i}(2,2:end))';
    
    if ymin+sqlength > nxl
        IN = IN(1:nxl-ymin+1,:);
    end
    if xmin+sqlength > nil
        IN = IN(:,1:nil-xmin+1);
    end
    closuregrd(ymin:min(ymin+sqlength,nxl),xmin:min(xmin+sqlength,nil)) = closuregrd(ymin:min(ymin+sqlength,nxl),xmin:min(xmin+sqlength,nil))+IN;    
end
closuregrd = spones(closuregrd);
closuregrd = closuregrd(1:nxl,1:nil);
end