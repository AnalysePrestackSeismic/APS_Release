function [wave] = ormsby(t_wave,t_sample)

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

    %required ormsby input parameters
    ormsby_f1=input('Enter low-cut frequency: ');
    ormsby_f2=input('Enter low-pass frequency: ');
    ormsby_f3=input('Enter high-pass frequency: ');
    ormsby_f4=input('Enter high-cut frequency: ');

    %preallocate ormsby array
    ormsby_sym=(zeros(length(t_axis_sym),1));

    %preallocate sinc2 arrays
    sinc2_f1=(zeros(length(t_axis_sym),1));
    sinc2_f2=(zeros(length(t_axis_sym),1));
    sinc2_f3=(zeros(length(t_axis_sym),1));
    sinc2_f4=(zeros(length(t_axis_sym),1));

    %build sinc2 functions for ormsby wavelet
    for t=1:length(t_axis_sym)
        if t_axis_sym(t,1)==0 %avoids division by zero
            sinc2_f1(t,1)=1;
            sinc2_f2(t,1)=1;
            sinc2_f3(t,1)=1;
            sinc2_f4(t,1)=1;
        else
            sinc2_f1(t,1)=(sin(pi*ormsby_f1*t_axis_sym(t,1))/(pi*ormsby_f1*t_axis_sym(t,1)))*(sin(pi*ormsby_f1*t_axis_sym(t,1))/(pi*ormsby_f1*t_axis_sym(t,1)));        
            sinc2_f2(t,1)=(sin(pi*ormsby_f2*t_axis_sym(t,1))/(pi*ormsby_f2*t_axis_sym(t,1)))*(sin(pi*ormsby_f2*t_axis_sym(t,1))/(pi*ormsby_f2*t_axis_sym(t,1)));
            sinc2_f3(t,1)=(sin(pi*ormsby_f3*t_axis_sym(t,1))/(pi*ormsby_f3*t_axis_sym(t,1)))*(sin(pi*ormsby_f3*t_axis_sym(t,1))/(pi*ormsby_f3*t_axis_sym(t,1)));
            sinc2_f4(t,1)=(sin(pi*ormsby_f4*t_axis_sym(t,1))/(pi*ormsby_f4*t_axis_sym(t,1)))*(sin(pi*ormsby_f4*t_axis_sym(t,1))/(pi*ormsby_f4*t_axis_sym(t,1)));
        end
    end

    %build ormsby wavelet
    for t=1:length(t_axis_sym)
        ormsby_sym(t,1)=((((pi*(ormsby_f4^2))/(ormsby_f4-ormsby_f3))*sinc2_f4(t,1))-(((pi*(ormsby_f3^2))/(ormsby_f4-ormsby_f3))*sinc2_f3(t,1)))-((((pi*(ormsby_f2^2))/(ormsby_f2-ormsby_f1))*sinc2_f2(t,1))-(((pi*(ormsby_f1^2))/(ormsby_f2-ormsby_f1))*sinc2_f1(t,1)));
    end

    %normalise ormsby wavelet
    ormsby_sym=ormsby_sym/(max(ormsby_sym));

    %calculate ormsby fft
    ftdata_ormsby_sym=fft(ormsby_sym);
    ftdata_norm_ormsby_sym=abs(ftdata_ormsby_sym)/(max(abs(ftdata_ormsby_sym)));
    
    if ftdata_norm_ormsby_sym((round(w_ns/2)),1) > 0.0125
        fprintf('Warning: non-zero amplitude at Nyquist frequency. Reduce samlping interval or central frequency to avoid aliasing.\n')
        cont=input('Do you wish to continue? 0 = No, 1 = Yes: ');
    else
        cont=1;
    end
    
    if cont==0
        return
    end

    if ploton==1 
        %display ormsby wavelet
        figure(2);
        subplot(2,1,1);
        plot(t_axis_sym,ormsby_sym);
        xlabel('Time (s)');
        ylabel('Amplitude');
        title(sprintf('%d-%d-%d-%dHz Ormsby Wavelet Time Series',ormsby_f1,ormsby_f2,ormsby_f3,ormsby_f4));
        axis tight;
        %display amplitude spectrum of ricker wavelet
        subplot(2,1,2);
        plot(f_axis,ftdata_norm_ormsby_sym((1:(round(w_ns/2))),1));
        xlabel('Frequency (Hz)');
        ylabel('Amplitude');
        title('Normalised Wavelet Amplitude Spectrum');
        axis tight;
    end

    wave=[t_axis_sym ormsby_sym f_axis_long ftdata_ormsby_sym];
end