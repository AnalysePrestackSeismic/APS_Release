clear all
test = fread(fopen('flat_wheeler_seabed-deep-prop-8_vol_time.bin'),'float32');
test = reshape(test,466,[]);
test = test';
[C,h] = contour(reshape(test(:,1),645,[]),25);

loopcount = 1;
ifcount = 1;
jump = 1;
nclosures = 0;

while jump < length(C)
    Z(loopcount,1) = C(1,jump);
    if C(1,jump+1) == C(1,jump+C(2,jump)) && C(2,jump+1) == C(2,jump+C(2,jump))
        nclosures = nclosures+1;
        closures{ifcount,1} = C(:,jump:jump+C(2,jump));
        ifcount = ifcount+1;
    end
    jump = jump+C(2,jump)+1;
    loopcount = loopcount+1;
end

fclose all;

grd = zeros(645,295);

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
    
    grd(ymin:ymin+sqlength,xmin:xmin+sqlength) = IN;
end

    