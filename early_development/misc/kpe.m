function [phase] = kpe(data,nil_kpe,nxl_kpe,nt_kpe)
% Phase estimation using maximum kurtosis method.

dt = 0.004;
data = reshape(data{1},[],1);
phi = (-90:1:90)*(pi/180);
kurt = zeros(length(phi),1);

% Wavelet parameters
tw = 0.15;                  % wavelet duration
nw_h = round(tw/dt/2);      % half number of wavelet samples
nw = 2*nw_h-1;              % total number of wavelet samples

% Wiener damping constant
damp = 3;                   % damping constant for Wiener filtering

% Inject small amount of noise to stabilise algorithm
snr = 20;                                 % S/N ratio
noise = randn(length(data),1);            % Gaussian noise
noise = noise*std(data)/std(noise)/snr;   % rescale noise to have correct S/N
data_nsy = data+noise;                    % add noise

data_nsy = reshape(data_nsy,nt_kpe,(nxl_kpe*nil_kpe));

% Create initial zero phase wavelet from power spectrum
ftdata_av = mean(abs(fft(data_nsy)),2);                         % average amplitude spectra
w0 = ifft(ftdata_av/sqrt(nt_kpe),'symmetric');                      % zero phase estimate
w0 = (2/sqrt(pi))*w0(1:nw_h).*cos(.5*pi*(0:nw_h-1)'/(nw_h-.5)); % cos taper edge
sn0 = w0(1)+2*sum(w0(2:end).*(-1).^(1:nw_h-1)');
w0(1) = w0(1)-sn0;                                              % -> w0 below has a double roots -1
w0 = [w0(nw_h:-1:2); w0];                                       % construct rest

% Estimate noise level from data spectrum between f_nyquist/2 and f_nyquist
vn0 = mean(ftdata_av(round(nt_kpe/4)+1:round(nt_kpe/2)+1).^2/nt_kpe);

% Normalise initial wavelet
w0 = w0/sqrt(vn0); 

% Zero-phase Wiener filter input data to improve phase estimate
ftw = fft([w0; zeros(nt_kpe,1)]);                                           % use zero padding to prevent wrap around
ftdata = fft([data_nsy; zeros(nw,nil_kpe*nxl_kpe)]);                            % move into frequency domain
G_est_damped = conj(ftw)./(abs(ftw).^2 + damp*damp)/sqrt(vn0);          % build Wiener filter
r_est_damped = ifft(G_est_damped*ones(1,nil_kpe*nxl_kpe).*ftdata, 'symmetric');    % apply Wiener filter and move back to time domain
r_est_damped_hilbert=imag(hilbert(r_est_damped));                       % apply hilbert transform

data = reshape(r_est_damped,[],1);
hilb_data = reshape(r_est_damped_hilbert,[],1);
for l = 1:length(phi)
   data_rot = data*cos(phi(l))+hilb_data*sin(phi(l));
   kurt(l,1) = bliss_kurtosis(data_rot);
end
[~,index] = max(kurt);
phase = phi(index)*180/pi;
end

