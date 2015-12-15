for i=1:size(lines,2)
    lines(i).len=sqrt((lines(i).point1(1)-lines(i).point2(1))^2+...
        (lines(i).point1(2)-lines(i).point2(2))^2);
    lines(i).depth=lines(i).point2(2)-lines(i).point1(2);
    len(i)=round(lines(i).len);
    depth(i)=round(lines(i).depth);
end

figure, hist(depth)