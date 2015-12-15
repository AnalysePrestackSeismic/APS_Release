function [f,F] = ricker_wavelet_function(fp)

s_rate = 0.001; % sample rate

% ricker wavelet in time domain
m = 1;
for i=1:1:length(fp)
    for j=-0.2:s_rate:0.2
        f(m,i) = (1-2*pi^2*fp(i)^2*j^2)*exp(-pi^2*fp(i)^2*j^2);
        m = m + 1;
    end
    m = 1;
end

% work out approximate size of ricker wavelets
for i=1:1:length(fp)
    T(i) = (sqrt(6)/pi*fp(i));
end

% ricker wavelet in frequency domain
m = 1;
n_max = 1/s_rate;
for i=1:1:length(fp)
    for n=0:1:n_max-1
        F(i,m) = (2/sqrt(pi))*(n^2/fp(i)^3)*exp(-n^2/fp(i)^2);
        m = m + 1;
    end
    m = 1;
end