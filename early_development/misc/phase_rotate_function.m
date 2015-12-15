function [rotated_time_series] = phase_rotate_function(time_series,phase)
% Phase rotation function

rotated_time_series = cos(-(phase*pi()/180))*time_series+sin(-(phase*pi()/180))*imag(hilbert(time_series));
end

