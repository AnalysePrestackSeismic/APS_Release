
clear all

% Local kurtosis

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

f = cos(-deg2rad(40))*f+sin(-deg2rad(40))*imag(hilbert(f));
s5 = filter2(f,r5,'same');

a = s1+s2+s3+s4+s5;
data = a;
i = 1;

N = 500;
tol = 0.00001;
lambda = 0.1;
% H = eye(length(a),length(a)); % should be same as w above;

c1 = ones(length(a));
%c2 = ones(length(d));
H = zeros(length(a),length(a));
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


% create seismic
for i_phase = -180:10:180;
    a = cos(-deg2rad(i_phase))*a+sin(-deg2rad(i_phase))*imag(hilbert(a));

    a = a.^2;
    A = diag(a);
    b = ones(length(a),1);
    B = diag(b);

    [c1] = conjugate_gradient(A,H,b,lambda,tol,N);

    [c2] = conjugate_gradient(B,H,a,lambda,tol,N);

    gamma(:,i) = 1./(c1.*c2).^2;
    i = i + 1;
end

figure
subplot(4,1,1); imagesc(data); axis tight
subplot(4,1,2); imagesc(b); axis tight
subplot(4,1,3); imagesc(c1); axis tight
subplot(4,1,4); imagesc(c2); axis tight
subplot(5,1,5); imagesc(gamma); axis tight
% subplot(4,1,4); imagesc(H'*H); axis tight