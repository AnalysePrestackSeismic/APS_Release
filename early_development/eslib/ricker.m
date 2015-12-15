function [wave] = ricker(t_wave,t_sample)

    ploton=1;

    %required global wavelet parameters
    w_length=t_wave;
    w_sample=t_sample;

    %calculate number of wavelet samples
    if floor(w_length/w_sample)==w_length/w_sample
        w_ns=w_length/w_sample;
    else
        w_ns=(floor(w_length/w_sample))+1;
    end

    %calculate sampling frequency
    f_sample = 1/w_sample;

    %build time axis
    %t_axis=(0:w_sample:w_length)';
    if floor(w_ns/2)==w_ns/2
        t_axis_sym=(-((w_ns/2-1)*w_sample):w_sample:((w_ns/2-1)*w_sample))';
    else
        t_axis_sym=(-((floor(w_ns/2))*w_sample):w_sample:((floor(w_ns/2))*w_sample))';
    end

    %build frquency axis
    if floor(w_ns/2)==w_ns/2
        f_axis=((0:1:((w_ns)/2)-1))*(f_sample/(w_ns-2))';
        f_axis_long=(0:1:w_ns-2)'*(f_sample/(w_ns-2));
    else
        f_axis=((0:1:(floor((w_ns)/2)))*(f_sample/(w_ns-1)))';
        f_axis_long=(0:1:w_ns-1)'*(f_sample/(w_ns-1));
    end

    %required ricker input parameters
    ricker_cf=input('Enter the central frequency: ');

    %preallocate ricker array
    ricker_sym=(zeros(length(t_axis_sym),1));

    %build ricker wavelet
    for t=1:length(t_axis_sym)
        ricker_sym(t,1)=(1-(2*(pi^2)*(ricker_cf^2)*(t_axis_sym(t,1)^2)))*exp(-(pi^2)*(ricker_cf^2)*(t_axis_sym(t,1)^2));
    end

    %calculate ricker fft
    ftdata_ricker_sym=fft(ricker_sym);
    ftdata_norm_ricker_sym=abs(ftdata_ricker_sym)/(max(abs(ftdata_ricker_sym)));
    
    if ftdata_norm_ricker_sym((round(w_ns/2)),1) > 0.0125
        fprintf('Warning: non-zero amplitude at Nyquist frequency. Reduce samlping interval or central frequency to avoid aliasing.\n')
        cont=input('Do you wish to continue? 0 = No, 1 = Yes: ');
    else
        cont=1;
    end
    
    if cont==0
        return
    end

    if ploton==1 
        %display ricker wavelet
        figure(1);
        subplot(2,1,1);
        plot (t_axis_sym,ricker_sym);
        xlabel('Time (s)');
        ylabel('Amplitude');
        title(sprintf('%dHz Ricker Wavelet Time Series',ricker_cf));
        axis tight;
        %display amplitude spectrum of ricker wavelet
        subplot(2,1,2);
        plot(f_axis,ftdata_norm_ricker_sym((1:(round(w_ns/2))),1));
        xlabel('Frequency (Hz)');
        ylabel('Amplitude');
        title('Normalised Wavelet Amplitude Spectrum');
        axis tight;
    end

    wave=[t_axis_sym ricker_sym f_axis_long ftdata_ricker_sym];

end