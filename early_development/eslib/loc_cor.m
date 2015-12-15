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

N = 1000;
tol = 0.0000000001;
lambda = 0;%0.00000001; %norm(a,2);
%O = ones(length(a)); %,length(a)); % should be same as w above;
%J = zeros(length(a),length(a));
%S = spdiags(O,0,J);

v = 1; %50:-1:0; %10:-1:0; %2:-1:0; %30:-1:0;
%v = v./length(v);
c3 = zeros(length(a),1);
c3(1:length(v)) = v;
%c1 = c1./length(v);
%c2 = ones(length(d));
%S = zeros(length(a),length(a));
%v = [-10:1:10];

%S = spdiags(c1,v,S);
%S = S./length(v);
%S = spdiags(c1,0,S);
S = toeplitz(c3);
H = S.^(1/2);
%H = spdiags(c1,[-1 1],H);
%H = spdiags(c1,[-2 2],H);
%H = spdiags(c1,[-3 3],H);
%H = spdiags(c1,[-4 4],H);
%H = spdiags(c1,[-5 5],H);
%H = spdiags(c1,[-6 6],H);
%H = spdiags(c1,[-7 7],H);
%H = spdiags(c1,[-8 8],H);
%H = spdiags(c1,[-9 9],H);


% create seismic
i = 1;
for i_phase = -180:5:180;
    f = cos(-(i_phase*pi()/180))*a+sin(-(i_phase*pi()/180))*imag(hilbert(a));
    %j(:,i) = a;
    hx = hilbert(f);
    % calculating the INSTANTANEOUS AMPLITUDE (ENVELOPE)
    b = abs(hx);
    %b = a;
    B = diag(b);
    % b = -a;
    %b = ones(length(a),1);
    %B = diag(b);

    d = f.^2;
    D = diag(f);
     
    [c1] = conjugate_gradient(D,H,b,lambda,tol,N);

    [c2] = conjugate_gradient(B,H,d,lambda,tol,N);

    gamma(:,i) = c1.*c2; % 1./(c1.*c2).^2;
    %gamma2(:,i) = 1./(gamma(:,i).^2);
    i = i + 1;
end 
figure(1)
subplot(5,1,1); plot(a); hold on; p = plot(b); set(p,'Color','red'); hold off; axis tight
subplot(5,1,2); plot(b); axis tight
subplot(5,1,3); plot(c1); axis tight
subplot(5,1,4); plot(c2); axis tight
subplot(5,1,5); plot(gamma); axis tight
% subplot(4,1,4); imagesc(H'*H); axis tight

figure(2)
imagesc(S);

figure(3)
%clims = [-0.05 0.05];
imagesc(gamma);

figure(5);
subplot(2,1,1); plot(gamma(450,:));
%subplot(2,1,2); plot(gamma2(450,:));