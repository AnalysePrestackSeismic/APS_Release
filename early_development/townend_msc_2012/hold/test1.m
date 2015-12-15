%RGB = imread('gantrycrane.png');
%I  = rgb2gray(RGB); % convert to intensity
%BW = edge(I,'canny'); % extract edges
[H,T,R] = hough_(traces.data,'RhoResolution',0.5,'Theta',-90:0.5:89.5);

% display the original image
subplot(2,1,1);
imshow(RGB);
title('gantrycrane.png');

% display the hough matrix
subplot(2,1,2);
imshow(imadjust(mat2gray(H)),'XData',T,'YData',R,...
    'InitialMagnification','fit');
title('Hough transform of gantrycrane.png');
xlabel('\theta'), ylabel('\rho');
axis on, axis normal, hold on;
colormap(hot);