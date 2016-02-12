function [peaks_out, peaks_idx ] = find_peaks( input_vector,threshold,gap )
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

%FIND_PEAKS reduce a vector to a sparse set of peaks
%   Detailed explanation goes here

% input should be a vector

if min(size(input_vector)) ~= 1 || ndims(input_vector) ~= 2
    error('Input must be vector.')
end
    

peaks_out = zeros(size(input_vector));

while max(input_vector)>threshold % threshold is a user parameter
    [~,idx] = max(input_vector);
    peaks_out(idx) = input_vector(idx);
    start_mask=max(idx-gap,1);
    end_mask=min(idx+gap,length(peaks_out));
    input_vector(start_mask:end_mask) = 0;
end

peaks_idx = find(peaks_out);


end

