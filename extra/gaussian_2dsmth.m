function [ smth_grid ] = gaussian_2dsmth( input_grid,hsmth,vsmth )
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
%
%   Create and apply a 2d gaussian smoother to a grid
%   
%   input_grid : grid to be smoothed
%   hmsth : size of horizontal smoother (diameter)
%   vsmth : size of vertical smoother (diameter)
%
%   input is padded using nearest neighbour prior to smoothing
%   with an elliptical 2D gaussian function 
%
%   DanB 22/01/2015

% force smoother lengths to be odd

hsmth = hsmth+(1-mod(hsmth,2));
vsmth = vsmth+(1-mod(vsmth,2));

    

[R C] = ndgrid(1:vsmth,1:hsmth);

sigma = max(hsmth,vsmth);

filter = gauss_calc(C,R,sigma,hsmth,vsmth);
filter = filter./(sum(sum(filter))); % normalise the filter

hpad = (hsmth-1)/2; vpad = (vsmth-1)/2;

% crappy coding of padding, there must be a more elegant way

padgrid1 = ones(vpad,hpad).*input_grid(1,1);
padgrid2 = repmat(input_grid(1,:),vpad,1);
padgrid3 = ones(vpad,hpad).*input_grid(1,end);
padgrid4 = repmat(input_grid(:,1),1,hpad);
padgrid5 = repmat(input_grid(:,end),1,hpad);
padgrid6 = ones(vpad,hpad).*input_grid(end,1);
padgrid7 = repmat(input_grid(end,:),vpad,1);
padgrid8 = ones(vpad,hpad).*input_grid(end,end);

padgrid = vertcat([padgrid1,padgrid2,padgrid3],[padgrid4,input_grid,padgrid5],[padgrid6,padgrid7,padgrid8]);

smth_grid=convn(padgrid,filter,'valid');

end

function val = gauss_calc(x,y,sigma,x_max,y_max)

y = y*(x_max/y_max);

x0 = (x_max+1)/2; 
y0 = ((y_max+1)/2)*(x_max/y_max);

exponent = ((x-x0).^2 + (y-y0).^2)./(2*sigma^2);
amplitude = 1 / (2 * sqrt(2*pi));

val = amplitude  * exp(-exponent);

end

