
function [frequency, maximum, gradient, index] = spec_gathers(gather, freq_range, time_range)

% function to pick peaks of spectral decomposition gathers
% Input: 
%       - a pre-stack spectral decomposition gather
%       - column is frequency
% Output:
%       -
%       - 
%       - change in frequency with time

% NOTE - might need to add time step

[samples,frequencies] = size(gather);

[maximum, index] = max(gather,[],2);

frequency = freq_range(index)';

for i=1:1:samples-1
    gradient(i) = (frequency(i+1)-frequency(i))/time_range(3);
end

gradient = gradient';