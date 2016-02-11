function c = curve_fitting(x,y,n)
%% ------------------ Disclaimer  ------------------
% 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) makes no representation or warranty, express or implied, in 
% respect to the quality, accuracy or usefulness of this repository. The code
% is this repository is supplied with the explicit understanding and 
% agreement of recipient that any action taken or expenditure made by 
% recipient based on its examination, evaluation, interpretation or use is 
% at its own risk and responsibility.
% 
% No representation or warranty, express or implied, is or will be made in 
% relation to the accuracy or completeness of the information in this 
% repository and no responsibility or liability is or will be accepted by 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) in relation to it.
%% ------------------ License  ------------------ 
% GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
%% github
% https://github.com/AnalysePrestackSeismic/
%% ------------------ FUNCTION DEFINITION ---------------------------------
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