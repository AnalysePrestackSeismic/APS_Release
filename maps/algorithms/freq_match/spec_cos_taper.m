function [ freq_pairs ] = spec_cos_taper( f1,f2,f3,f4 )
%generate freq, amplitude pairs from 4 corner points
%using supplied 4 corner frequencies
%note that this is hardwired for 1 Hz intervals which may not be good
%enough for low frequencies

freq_pairs = zeros(f4,2); % initialise output

freq_pairs(:,1) = [1:f4]; % set frequency values to range of corner points

% frequencies between f1 and f2 set to cosine taper up from 0 to 1
% frequencies between f2 and f3 are set to 1
% frequencies between f3 and f4 set to cosine taper down from 1 to 0

freq_pairs(f1:f2,2) = fliplr(0.5*(1+cos(-pi/(f2-f1).*([f1:f2]-f1))));
freq_pairs(f2:f3,2) = 1;
freq_pairs(f3:f4,2) = 0.5*(1+cos(-pi/(f4-f3).*([f3:f4]-f3)));

% plot the output

figure; plot(freq_pairs(:,1),freq_pairs(:,2))


end

