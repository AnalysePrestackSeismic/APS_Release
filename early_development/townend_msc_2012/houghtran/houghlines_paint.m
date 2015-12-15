for i=1:size(lines,2)
    lines(i).len=sqrt((lines(i).point1(1)-lines(i).point2(1))^2+...
        (lines(i).point1(2)-lines(i).point2(2))^2);
    lines(i).dep=lines(i).point2(2)-lines(i).point1(2);
end

[tmp ind]=sort([lines.dep]);
lines=lines(fliplr(ind));

wid=[2,2];
Hslice=zeros(size(traces.data(1:500,:)));
for i=1:size(lines,2);
    for j=0:(1/lines(i).len):1
        loc=round(lines(i).point1+j*(lines(i).point2-lines(i).point1));
        ind1=loc(1)-wid(1):loc(1)+wid(1);
        ind2=loc(2)-wid(2):loc(2)+wid(2);
        Hslice(ind1,ind2)=lines(i).dep*traces.data(ind1,ind2)/10;
    end
end