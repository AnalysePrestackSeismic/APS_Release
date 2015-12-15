clear all
% Local instaneous frequency

% make time series
t = 0:0.002:1;
d = sin(2*pi*25*t.*exp(t));
data = d;

% 
N = 500;
tol = 0.00001;
lambda = 1;

dif_d = diff(d); % differentiate d
hilb_d = imag(hilbert(d)); % take the hilbert transform of d (seems no point taking real part since it is just d)?

L = diag(d(1:500).^2+hilb_d(1:500).^2); % denominator f(t)^2+hilb_d^2

dif_hilb_d = diff(hilb_d); % differentiate hilbert of d
d = (d(1:500).*dif_hilb_d)-(dif_d.*hilb_d(1:500)); % numerator f(t)h'(t)+f'(t)h(t)
d = d';
w = inv(L)*d;

% make shaping matrix
H = eye(length(d),length(d)); % should be same as w above;

%c1 = ones(length(d));
%c2 = ones(length(d));
%H = zeros(length(d),length(d));
%v = [-10:1:10];
% H = spdiags(c1,0,H);
% H = spdiags(c1,[-1 1],H);
% H = spdiags(c1,[-2 2],H);
% H = spdiags(c1,[-3 3],H);
% H = spdiags(c1,[-4 4],H);
% H = spdiags(c1,[-5 5],H);
% H = spdiags(c1,[-6 6],H);
% H = spdiags(c1,[-7 7],H);
% H = spdiags(c1,[-8 8],H);
% H = spdiags(c1,[-9 9],H);

[m] = conjugate_gradient(L,H,d,lambda,tol,N);

figure
subplot(4,1,1); plot(data); axis tight
subplot(4,1,2); plot(w); axis tight
subplot(4,1,3); plot(m); axis tight