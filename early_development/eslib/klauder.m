function [wave] = klauder(t_wave,t_sample)

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

    %required klauder input parameters
    klauder_f1=input('Enter low-cut frequency: ');
    klauder_f2=input('Enter high-cut frequency: ');
    klauder_t=input('Enter signal duration: ');
    small_no=0.0000000000000001; %avoids division by zero

    %calculate some intermediate parameters
    k=(klauder_f2-klauder_f1)/klauder_t;
    f0=(klauder_f2+klauder_f1)/2;

    %preallocate klauder array
    klauder_sym=zeros(length(t_axis_sym),1);

    %build klauder wavelet
    for t=1:length(t_axis_sym)
        if t_axis_sym(t,1)==0
            klauder_sym(t,1)=(sin(pi*k*small_no*(klauder_t-small_no)))/((pi*k*small_no)*exp(2*pi*i*f0*small_no));
        else
            klauder_sym(t,1)=(sin(pi*k*t_axis_sym(t,1)*(klauder_t-t_axis_sym(t,1))))/((pi*k*t_axis_sym(t,1))*exp(2*pi*i*f0*t_axis_sym(t,1)));
        end
    end

    %calculate klauder fft
    ftdata_klauder_sym=fft(real(klauder_sym));
    ftdata_norm_klauder_sym=abs(ftdata_klauder_sym)/(max(abs(ftdata_klauder_sym)));
    
    if ftdata_norm_klauder_sym((round(w_ns/2)),1) > 0.0125
        fprintf('Warning: non-zero amplitude at Nyquist frequency. Reduce samlping interval or central frequency to avoid aliasing.\n')
        cont=input('Do you wish to continue? 0 = No, 1 = Yes: ');
    else
        cont=1;
    end
    
    if cont==0
        return
    end
 

    if ploton==1 
        %display klauder wavelet
        figure(3);
        subplot(2,1,1);
        plot(t_axis_sym,(real(klauder_sym)/max(real(klauder_sym))));
        xlabel('Time (s)');
        ylabel('Amplitude');
        title(sprintf('%d-%dHz Klauder Wavelet Time Series (T=%ds)',klauder_f1,klauder_f2,klauder_t));
        axis tight;
        %display amplitude spectrum of ricker wavelet
        subplot(2,1,2);
        plot(f_axis,ftdata_norm_klauder_sym((1:(round(w_ns/2))),1));
        xlabel('Frequency (Hz)');
        ylabel('Amplitude');
        title('Normalised Wavelet Amplitude Spectrum');
    end

    wave=[t_axis_sym real(klauder_sym)/max(real(klauder_sym)) f_axis_long ftdata_klauder_sym];
end
