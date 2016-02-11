function [ smth_grid ] = gaussian_1dsmth( input_grid,smth_size )
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
%   Create and apply a 1-d gaussian smoother to a grid

t = linspace(-1,1,smth_size)';
filttraces = zeros(length(t),1);
a = 1;
filttraces = sqrt(pi)/a*exp(-(pi*t/a).^2);
filttraces = filttraces/sum(filttraces);
filttraces = filttraces';

padlen = floor(size(filttraces,2)/2);
input_grid=[repmat(input_grid(:,1),1,padlen) input_grid repmat(input_grid(:,size(input_grid,2)),1,padlen)];
smth_grid=convn(input_grid,filttraces,'valid');



end

