img=squeeze(lines);
rotI=img;
BW = rotI;
[H,T,R] = hough(img);
figure('Color',[1 1 1]); subplot(2,6,[1,2,3,4,7,8,9,10]);
imagesc(H,'XData',T,'YData',R);%,'InitialMagnification','fit');
c=colorbar;
ylabel(c,'Accumulations');
xlabel('\theta (degrees)'), ylabel('\rho (pixels)');
axis on, axis normal, hold on;
P  = houghpeaks(H,3,'threshold',ceil(0.05*max(H(:))),'nhoodsize',2.*round(((size(H)/10)+1)/2)-1);
x = T(P(:,2));
y = R(P(:,1));
plot(x,y,'s','color','red');

% Find lines and plot them
lines = houghlines(img,T,R,P,'FillGap',4,'MinLength',5);
subplot(2,6,[5,6,11,12]); imagesc(BW), hold on
daspect([1 1 1])
max_len = 0;
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:,1),xy(:,2),'LineWidth',1,'Color','green');
    
    % plot beginnings and ends of lines
    plot(xy(1,1),xy(1,2),'x','LineWidth',1,'Color','yellow');
    plot(xy(2,1),xy(2,2),'x','LineWidth',1,'Color','red');
    
    % determine the endpoints of the longest line segment
    len = norm(lines(k).point1 - lines(k).point2);
    if ( len > max_len)
        max_len = len;
        xy_long = xy;
    end
end

% highlight the longest line segment
%plot(xy_long(:,1),xy_long(:,2),'LineWidth',1,'Color','black');