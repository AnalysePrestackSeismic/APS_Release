clear all

sparsity = 0.9;
skew = 0.4;
ai(1,1) = 1500;
ns = 4000;
sr = 0.004;


for ii = 1:ns
    if ii == 100
        r(ii,1) = 4;
    elseif ii > 100
        if rand>sparsity
            if rand>skew
                r(ii,1) = abs(randn/(ii^0.1));
            else
                r(ii,1) = -abs(randn/(ii^0.1));
            end
        else
            r(ii,1) = 0;
        end
    end
end

r = 0.25*r/max(abs(r));

for ii = 1:ns-1
    ai(ii+1,1) = ai(ii,1)*(r(ii,1)+1)/(1-r(ii,1));
end

w = ricker(35,sr);
nsw = length(w);
hwl = floor(nsw/2);

s = conv(r,w);

noise = randn(ns+nsw-1,1);
noise = noise/norm(noise);
snr = 10;
sn = s + (norm(s)*noise/norm(noise))/snr;

absftw = abs(fft(w,ns+nsw-1));
pw = 0*max(absftw)/1e1;

absftwpw = zeros(ns+nsw-1,1);
absftwpw(absftw>pw) = absftw(absftw>pw).^2;
absftwpw(absftw<=pw) = pw.^2;

rinv = ifft((fft(conj(w),ns+nsw-1).*fft(s,ns+nsw-1))./(fft(w,ns+nsw-1).*fft(w,ns+nsw-1)),'symmetric');
rinv = rinv(1:end-nsw,1);

% ghost

wg = conv(w,[zeros(nsw+6,1);-0.9;zeros(nsw-6,1)],'same');
wg = w+wg;

sg = conv(r,wg);

sng = sg + (norm(sg)*noise/norm(noise))/snr;

absftwg = abs(fft(wg,ns+nsw-1));
pwg = 0*max(absftwg)/1e1;

absftwpwg = zeros(ns+nsw-1,1);
absftwpwg(absftwg>pwg) = absftwg(absftwg>pwg).^2;
absftwpwg(absftwg<=pwg) = pwg.^2;

rinvg = ifft((fft(conj(wg),ns+nsw-1).*fft(sg,ns+nsw-1))./(fft(wg,ns+nsw-1).*fft(wg,ns+nsw-1)),'symmetric');

rinvg = rinvg(1:end-nsw,1);

f = 0:250/4033:250;

figure(1);
subplot(3,1,1)
plot(r)
axis([0 4000 -0.1 0.15])
subplot(3,1,2)
plot(rinv-0.0015)
axis([0 4000 -0.1 0.15])
subplot(3,1,3)
plot(rinvg)
axis([0 4000 -0.1 0.15])

% decon operators
decon = fft(conj(w),ns+nsw-1)./(fft(w,ns+nsw-1).*fft(w,ns+nsw-1));
decong = fft(conj(wg),ns+nsw-1)./(fft(wg,ns+nsw-1).*fft(wg,ns+nsw-1));
deconpg = unwrap(angle(decong))-(0.068*2*pi*f')-2*pi;
pg = unwrap(angle(fft(wg,4034)))+(0.068*2*pi*f');

figure(2);
subplot(3,2,1)
plot((-17:1:17)*4,w)
axis([-68 68 -1 1.5])
subplot(3,2,2)
plot((-17:1:17)*4,wg)
axis([-68 68 -1 1.5])
subplot(3,2,3)
plot(f,abs(fft(w,4034)))
hold all
plot(f,zeros(4034,1))
axis([0 125 -1.5 5])
subplot(3,2,4)
plot(f,abs(fft(wg,4034)))
hold all
plot(f,pg)
axis([0 125 -1.5 5])
subplot(3,2,5)
plot(f,abs(decon))
hold all
plot(f,zeros(4034,1))
axis([0 125 -1.5 5])
subplot(3,2,6)
plot(f,abs(decong))
hold all
plot(f,deconpg)
axis([0 125 -1.5 5])