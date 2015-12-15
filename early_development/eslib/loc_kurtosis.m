% Phase estimation based on correlation of envelope and signal amplitude
% squared
% Reference: 2010. Fomel, S. and van der Baan, M. Local similarity with envelope. SEG Denver 2010 Annual Meeting 
% 
% Authors: James Selvage and Jonathan Edgar (2011)
clear all

% generate a test ricker wavelet
[f] = ricker_wavelet_function(25);

% number of samples
n_samp = 501;

% this controls the spike location, magnitude and phase
spike_1 = [250 0.01 40];
spike_2 = [150 0.01 50];
spike_3 = [275 -0.01 90];
spike_4 = [300 0.01 0];
spike_5 = [450 0.01 -10];

% this makes the synthetic section
r1 = zeros(n_samp,1);

r1(spike_1(1)) = spike_1(2);
r2 = zeros(n_samp,1);
r2(spike_2(1)) = spike_2(2);
r3 = zeros(n_samp,1);
r3(spike_3(1)) = spike_3(2);
r4 = zeros(n_samp,1);
r4(spike_4(1)) = spike_4(2);
r5 = zeros(n_samp,1);
r5(spike_5(1)) = spike_5(2);

% create some seismic
f = cos(-(spike_1(3)*pi()/180))*f+sin(-(spike_1(3)*pi()/180))*imag(hilbert(f));
s1 = filter2(f,r1,'same');

f = cos(-(spike_2(3)*pi()/180))*f+sin(-(spike_2(3)*pi()/180))*imag(hilbert(f));
s2 = filter2(f,r2,'same');

f = cos(-(spike_3(3)*pi()/180))*f+sin(-(spike_3(3)*pi()/180))*imag(hilbert(f));
s3 = filter2(f,r3,'same');

f = cos(-(spike_4(3)*pi()/180))*f+sin(-(spike_4(3)*pi()/180))*imag(hilbert(f));
s4 = filter2(f,r4,'same');
% 
f = cos(-(-spike_5(3)*pi()/180))*f+sin(-(-spike_5(3)*pi()/180))*imag(hilbert(f));
s5 = filter2(f,r5,'same');

a = s1+s2+s3+s4+s5;

% a = a + 0.001*randn(length(a),1);

N = 1000; % maximum number of iterations
tol = 0.0000000001; % tolerance of inversion
lambda = 0; %0.00000001; %norm(a,2);

% local smoother to apply method over
% currently creates a triangle smoother
v = 1; %50:-1:0; %10:-1:0; %2:-1:0; %30:-1:0;
s_col = zeros(length(a),1);
s_col(1:length(v)) = v;
S = toeplitz(s_col);
H = S.^(1/2);

% create seismic
i = 1;
for i_phase = -180:5:180;
    s = cos(-(i_phase*pi()/180))*a+sin(-(i_phase*pi()/180))*imag(hilbert(a));
    
    % Testing correlation with similarity
    hx = hilbert(s);
    b1 = abs(hx);
    B1 = diag(b1);
    d1 = s.^2; % similarity
    D1 = diag(s);
    % D1 = diag(d1);
     
    % Testing correlation with a constant (kurtosis)
    b2 = ones(length(a),1);
    B2 = diag(b2);
    d2 = s.^2;
    D2 = diag(s);
    %D2 = diag(d2); % this I change
    
    % similarity
    [sim_c1] = conjugate_gradient(D1,H,b1,lambda,tol,N);
    [sim_c2] = conjugate_gradient(B1,H,d1,lambda,tol,N);
    
    sim_gamma(:,i) = (sim_c1.*sim_c2);
    
    % kurtosis
    [kur_c1] = conjugate_gradient(D2,H,b2,lambda,tol,N);
    [kur_c2] = conjugate_gradient(B2,H,d2,lambda,tol,N);
    
    kur_gamma(:,i) = (kur_c1.*kur_c2);

    i = i + 1;
end

figure(1)
subplot(4,1,1); plot(a); hold on; p = plot(b1); set(p,'Color','red'); hold off; axis tight  % p = plot(b2); set(p,'Color','blue'); hold off; axis tight
subplot(4,1,2); plot(sim_gamma(450,:)); hold on; j = plot(kur_gamma(450,:)); set(j,'Color','red'); hold off;  axis tight
subplot(4,1,3); imagesc(sim_gamma); axis tight
subplot(4,1,4); imagesc(kur_gamma); axis tight