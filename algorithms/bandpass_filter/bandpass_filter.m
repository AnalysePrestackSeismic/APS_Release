function  [o] =  bandpass_filter(d,dt,f1,f2,f3,f4);
%BP_Filter: Apply a band-pass filter to a group of traces.
%
%  [o] = bp_filter(d,dt,f1,f2,f3,f4);
%
%  IN   d:    data (columns are traces)
%       dt:   sampling interval in sec
%       f1:   freq. in Hz
%       f2:   freq. in Hz
%       f3:   freq. in Hz
%       f4:   freq. in Hz
%
%   ^
%   |     ___________
%   |    /           \   Amplitude spectrum
%   |   /             \
%   |  /               \
%   |------------------------>
%      f1 f2        f3 f4
%
%  OUT  o:    output  (columns are traces)
%
%  Example: 
%
%    d=linear_events; 
%    dout = bp_filter(d,0.004,1,3,30,40); 
%    wigb([d,dout]);
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
%  Author: M.D.Sacchi
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%
% modifications made to the tapers in freq and made to vectorise with
% bsxfun


 [nt,nx] = size(d);                                                     % Calculate the size of the data
 k = nextpow2(nt);                                                      % Calculate the smallest power of two greater than or equal to the total number of samples ina trace                         
 nf = 4*(2^k);                                                          % Design Length of the filter

 i1 = floor(nf*f1*dt)+1;                                                % Index of f1 in filter
 i2 = floor(nf*f2*dt)+1;                                                % Index of f2 in filter
 i3 = floor(nf*f3*dt)+1;                                                % Index of f3 in filter
 i4 = floor(nf*f4*dt)+1;                                                % Index of f4 in filter

 %up =  (1:1:(i2-i1))/(i2-i1);                                           % Design the filter taper ramp up
 up = ((sin(linspace((-(pi*0.87)/2),((pi)/2),(i2-i1))')+1)/2)';
 %down = (i4-i3:-1:1)/(i4-i3);                                           % Design the filter taper ramp down
 down = ((sin(linspace(((pi)/2),(-(pi*0.87)/2),(i4-i3))')+1)/2)';
 aux = [zeros(1,i1), up, ones(1,i3-i2), down, zeros(1,nf/2+1-i4) ];     % Design the whole filter joining the bits
 aux2 = fliplr(aux(1,2:nf/2));                                          % Flip the  filter

 c = 0;                                                                 % zero phase (could apply rotations as well)
 F = ([aux,aux2]');
 Phase = (pi/180.)*[0.,-c*ones(1,nf/2-1),0.,c*ones(1,nf/2-1)];          % Phase spectrum
 Transfer = F.*exp(-1i*Phase');


 D = fft(d,nf,1);

%  for k = 1:nx
%   Do(:,k) = Transfer.*D(:,k);
%  end

 Do(:,:) = bsxfun(@times,D(:,:),Transfer);
 
 
 o = ifft(Do,nf,1);

 o = real(o(1:nt,:));

 
