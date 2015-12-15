% create a synthetic seismogram for testing MCMC code

clear all

% user controls
fp = [50 45 30 80 25]; % peak frequency of ricker wavelets
n_samp = 800; % number of samples
s_rate = 0.004; % sample rate
sig = 5; % this controls the smoothing in the t-f panel
% spike locations and magnitude
spike_1 = [250 -0.2];
spike_2 = [300 0.3];
spike_3 = [400 -0.5];
spike_4 = [425 0.6];
spike_5 = [740 0.5];

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

refl = r1+r2+r3+r4+r5;

% create some seismic
s1 = filter2(f(:,1),r1,'same');
s2 = filter2(f(:,2),r2,'same');
s3 = filter2(f(:,3),r3,'same');
s4 = filter2(f(:,4),r4,'same');
s5 = filter2(f(:,5),r5,'same');

% create seismic
seis = s1+s2+s3+s4+s5;

% create time-frequency display
time_freq = zeros(n_samp,n_max);

for i=(spike_1(1)-round(T(1)/2)):1:(spike_1(1)+round(T(1)/2))
    time_freq(i,:) = time_freq(i,:) + F(1,:)*exp(-(spike_1(1)-i)^2/sig);
end

for i=(spike_2(1)-round(T(2)/2)):1:(spike_2(1)+round(T(2)/2))
    time_freq(i,:) = time_freq(i,:) + F(2,:)*exp(-(spike_2(1)-i)^2/sig);
end

for i=(spike_3(1)-round(T(3)/2)):1:(spike_3(1)+round(T(3)/2))
    time_freq(i,:) = time_freq(i,:) + F(3,:)*exp(-(spike_3(1)-i)^2/sig);
end

for i=(spike_4(1)-round(T(4)/2)):1:(spike_4(1)+round(T(4)/2))
    time_freq(i,:) = time_freq(i,:) + F(4,:)*exp(-(spike_4(1)-i)^2/sig);
end

for i=(spike_5(1)-round(T(5)/2)):1:(spike_5(1)+round(T(5)/2))
    time_freq(i,:) = time_freq(i,:) + F(5,:)*exp(-(spike_5(1)-i)^2/sig);
end

figure(1)
subplot(1,3,1);
plot(refl,(1:n_samp)*s_rate); set(gca,'YDir','reverse') 
ylabel('time (s)','FontSize',12)
xlabel('amplitude','FontSize',12)
axis tight;
xlim([-1 1]);

subplot(1,3,2);
plot(seis,(1:n_samp)*s_rate); set(gca,'YDir','reverse') 
xlabel('amplitude','FontSize',12)
set(gca,'YTick',[])
axis tight;
xlim([-1 1]);

subplot(1,3,3);
imagesc(time_freq);
xlabel('frequency (hz)','FontSize',12)
set(gca,'YTick',[])
title('time-frequency spectrum','FontSize',12)
axis tight;

% zero out low amplitudes
for i=1:1:length(seis)
   if abs(seis(i)) < 1e-4
       seis(i) = 0;
   end
end

% dlmwrite('decomp.dat',[100 100 seis'],'delimiter','\t','precision', 8)

% make time_freq take account of wavelet length
