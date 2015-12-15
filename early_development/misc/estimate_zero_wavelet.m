function wavelet = estimate_zero_wavelet(data,nt)
    wrms=1;

    Sft = mean(abs(fft(data)),2);

    wraw = ifft(Sft/sqrt(nt),'symmetric');

    wtmp = (2/sqrt(pi))*wraw(1:51).*cos(.5*pi*(0:50)'/(50.5));
    sn0 = wtmp(1)+2*sum(wtmp(2:end).*(-1).^(1:50)');
    wtmp(1) = wtmp(1)-sn0;
    wtmp = [wtmp(51:-1:2); wtmp];
    wtmp = wtmp/(norm(wtmp/wrms)/sqrt(length(wtmp)));
    if isnan(wtmp)
        wavelet = 0;
    else
        wavelet = wtmp;
    end
end   