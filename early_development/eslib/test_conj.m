clear all

t = 0:0.002:1;
d = sin(2*pi*50*t)+cos(2*pi*25*t);
d = d';

N = 10;
H = zeros(length(d));
H(1) = 3;
H(2) = 2;
H(3) = 1;
H = H./H(1);
H = H';
lambda = 0.1;
tol = 0.01;

L = diag(d);

[m] = conjugate_gradient(L,H,d,lambda,tol,N);