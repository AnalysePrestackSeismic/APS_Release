function [ smth_grid ] = gaussian_1dsmth( input_grid,smth_size )
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

