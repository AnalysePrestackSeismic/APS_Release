function c = curve_fitting(x,y,n)
%% Function to fit an nthe order polynomial for arrays x and y
% Use n=1 for straight line
% Function displays scatter of raw points and the fitted polynomial
% Function also displays the calcuate equation
%-----------------------------------------------------------
%Note: Always test for different values of n especially if problem is
%underdetermined or you are using this for extrapolation rather than
%inerpolation
%-------------------------------------------------------------
close all;
c=polyfit(x,y,n);

figure(1);
%Plot the result---------------------
xx=-1*ceil(0.5*max(x)):floor((max(x)-min(x))/50):ceil(1.5*max(x));
yy=polyval(c,xx);
scatter(x,y,20,'r');
hold on;
plot(xx,yy);
hold off;
%Print the polynomial on screen-----------------
fprintf ('The Polynomial function is:\t');
for i = 1:n
    fprintf('(%g)*x^%d+',c(i),n-i+1);
end

fprintf('(%g)\n',c(n+1));
end