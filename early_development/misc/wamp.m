function [wout] = wamp(data,nil_wamp,nxl_wamp,nt_wamp)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% Wavelet parameters
dt = 0.004;
tw = 0.15;                  % wavelet duration
nw_h = round(tw/dt/2);      % half number of wavelet samples

for l=1:length(data)

    % Create initial zero phase wavelet from power spectrum
    data_tmp = reshape(data{l,1},nt_wamp,nil_wamp*nxl_wamp);

    ftdata_av = mean(abs(fft(data_tmp)),2);                             % average amplitude spectra
    w = ifft(ftdata_av/sqrt(nt_wamp),'symmetric');                  % zero phase estimate
    w = (2/sqrt(pi))*w(1:nw_h).*cos(.5*pi*(0:nw_h-1)'/(nw_h-.5));   % cos taper edge
    sn0 = w(1)+2*sum(w(2:end).*(-1).^(1:nw_h-1)');
    w(1) = w(1)-sn0;                                                % -> w0 below has a double roots -1
    w = [w(nw_h:-1:2); w];                                          % construct rest
    wout(:,l) = w;
end

