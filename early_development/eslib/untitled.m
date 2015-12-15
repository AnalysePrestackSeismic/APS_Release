%t = 0:0.002:1;
%d = sin(2*pi*25*t.*exp(t));
[f] = ricker_wavelet_function(50);
n_samp = 501;
spike_1 = [50 -0.2];
spike_2 = [150 0.3];
spike_3 = [275 -0.5];
spike_4 = [350 0.6];
spike_5 = [450 0.5];

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
s1 = filter2(f,r1,'same');
s2 = filter2(f,r2,'same');
s3 = filter2(f,r3,'same');
s4 = filter2(f,r4,'same');
s5 = filter2(f,r5,'same');

% create seismic
d = s1+s2+s3+s4+s5;
%d = sin(2*pi*1*t);
%d2 = cos(2*pi*25*t);
%d = zeros(1,501);
%d(1,1:251) = d1(1:251);
%d(1,251:501) = d2(251:501);
data = d;

% d = d';
dif_d = diff(d);
hilb_d = imag(hilbert(d));
dif_hilb_d = diff(hilb_d);
L = diag(d(1:500).^2+hilb_d(1:500).^2);
d = (d(1:500).*dif_hilb_d)-(dif_d.*hilb_d(1:500));
%d = d';
w = inv(L)*d;

N = 1000;
%H = zeros(length(d),length(d));
%H = spdiags(ones(length(d));
%H = eye(length(d),length(d));
%H = ones(length(d));
c1 = ones(length(d));
%c2 = ones(length(d));
H = zeros(length(d),length(d));
%v = [-10:1:10];
H = spdiags(c1,0,H);
H = spdiags(c1,[-1 1],H);
H = spdiags(c1,[-2 2],H);
H = spdiags(c1,[-3 3],H);
H = spdiags(c1,[-4 4],H);
H = spdiags(c1,[-5 5],H);
H = spdiags(c1,[-6 6],H);
H = spdiags(c1,[-7 7],H);
H = spdiags(c1,[-8 8],H);
H = spdiags(c1,[-9 9],H);
%H = spdiags(c1,0,H);
%H = spdiags(c1.*2/3,[-1 1],H);
%H = spdiags(c1.*1/3,[-2 2],H);
%H = H./2;
% h1 = spdiags(c1,0,H);
% h2 = spdiags(c1,1,H);
% h3 = spdiags(c1,-1,H);
% H = h1+h2+h3;
%H = eye(length(d));

%H(1:N) = 1;
%H(1) = 3;
%H(2) = 2;
%H(3) = 1;
%H = diag(h);
%H(1) = 3;
%H(2) = 2;
%H(3) = 1;
%H = H./H(1);
%H = H';
lambda = 0.02;
tol = 0.001;

%L = diag(d);


figure
subplot(4,1,1); plot(data); axis tight
subplot(4,1,2); plot(w); axis tight
subplot(4,1,3); plot(m); axis tight
subplot(4,1,4); imagesc(H'*H); axis tight