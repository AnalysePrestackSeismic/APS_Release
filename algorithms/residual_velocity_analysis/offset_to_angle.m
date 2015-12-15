function [ offset_traces ] = offset_to_angle( job_meta,angles,vel_traces )
% calculates incidence angle for each offset and time
%   
% angles: array of angles to calculate offsets for
% vel_traces: 2D array of velocity traces
%
% returns offset_traces which is 2d array of angle vs sample
%

[nsamples ntraces] = size(vel_traces);

offset_traces = zeros(nsamples,size(offsets),ntraces);

angles_rad = 2*pi*angles/360;

% ======================================================================
%
% dix convert rms to interval
%
% replace with constrained inversion
%

vint(1,:) = vel_traces(1,:);
tt = (1:nsamples).*job_meta.srate;

for samp=2:nsamples
    
    vint(samp,:) = sqrt(((vel_traces(samp,:).^2*tt(samp)) - (vel_traces(samp-1,:).^2*tt(samp-1))) ./ (tt(samp)-tt(samp-1)));
    
end

% =====================================================================

% angle calculation

% sine theta = (offset * vint) / (t * vrms^2)

% => offset = (t * vrms^2 * sine theta) / vint


for ang_idx=1:size(angles)
    
    for trc_idx=1:ntraces
    
        offset_traces(:,ang_idx,trc_idx) = (tt(:) .* vel_traces(:,trc_idx) .* sin(angles_rad(ang_idx))) / vint(:,trc_idx);
    
    end
end


end

