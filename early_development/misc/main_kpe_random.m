%for dir=1:8
    
%cd (sprintf('%d',dir));

clear all
close all
clc

%Read data
file_in = input('''Input file name'' (segy) = ');                                        % input Segy file
endian = 'b';                                                                            % Segy endian
[Segy,SegyTraceHeaders,SegyHeader]=ReadSegyConstantTraceLength(file_in,'endian',endian); % calls SegyMAT to read Segy
%load (sprintf('%d.mat',dir));
%load synsect.mat synsect;
%Segy=synsect;


% Extract important header info
dt = SegyHeader.dt/1.e6;    % time sample size in Segy
%dt=0.004;
tmin = input('Input min time in seconds (enter 0 to start at beginning): ');
tmax = input('Input max time in seconds (enter 0 to finish at end): ');
[nt nrec] = size(Segy);
if tmin == 0
    smin=1;
else
    smin=round(tmin/dt);
end
if tmax == 0
    smax=nt;
else
    smax=round(tmax/dt);
end
if (tmax ~= 0) || (tmin ~=0)
    Segy=Segy(smin:smax,1:nrec);
end
[nt nrec] = size(Segy);     % number of time samples and number of records (traces)
fsample = 1./dt;            % sampling frequency

% Wavelet parameters
tw = 0.15;                   % wavelet duration
nw_h = round(tw/dt/2);      % half number of wavelet samples
nw = 2*nw_h-1;              % total number of wavelet samples

% Wiener damping constant
damp = 3;                   % damping constant for Wiener filtering

% Inject small amount of noise to stabilise algorithm
snr = 20;                                           % S/N ratio
for irec = 1:nrec
  noise = randn(nt,1);                              % Gaussian noise
  noise = noise*std(Segy(:,irec))/std(noise)/snr;   % rescale noise to have correct S/N
  sn = std(noise);                                  % noise standard deviation
  data_nsy(:,irec) = Segy(:,irec)+noise;            % add noise
end

% Create initial zero phase wavelet from power spectrum
ftdata_av = mean(abs(fft(data_nsy)),2);                         % average amplitude spectra
w0 = ifft(ftdata_av/sqrt(nt),'symmetric');                      % zero phase estimate
w0 = (2/sqrt(pi))*w0(1:nw_h).*cos(.5*pi*(0:nw_h-1)'/(nw_h-.5)); % cos taper edge
sn0 = w0(1)+2*sum(w0(2:end).*(-1).^(1:nw_h-1)');
w0(1) = w0(1)-sn0;                                              % -> w0 below has a double roots -1
w0 = [w0(nw_h:-1:2); w0];                                       % construct rest

% Estimate noise level from data spectrum between f_nyquist/2 and f_nyquist
vn0 = mean(ftdata_av(round(nt/4)+1:round(nt/2)+1).^2/nt);

% Normalise initial wavelet
w0 = w0/sqrt(vn0); 

% Zero-phase Wiener filter input data to improve phase estimate
ftw = fft([w0; zeros(nt,1)]);                                           % use zero padding to prevent wrap around
ftdata = fft([data_nsy; zeros(nw,nrec)]);                               % move into frequency domain
G_est_damped = conj(ftw)./(abs(ftw).^2 + damp*damp)/sqrt(vn0);          % build Wiener filter
r_est_damped = ifft(G_est_damped*ones(1,nrec).*ftdata, 'symmetric');    % apply Wiener filter and move back to time domain
r_est_damped_hilbert=imag(hilbert(r_est_damped));                       % apply hilbert transform

% Phase estimation
% Set parameters
nphase = 61;                  % number of phase estimates
dphase = 180/(nphase-1);      % phase sampling
sample_kurt=zeros(nphase,1);  % preallocate sample kurtosis array

% Loop parameters
iter_max = 400;             % maximum number of iterations (recommend >100)
ndata_min = 100000;         % minimum number of Segy samples to use (recommend >50000)
ndata_max = nt*nrec;        % maximum number of Segy samples to use (up to all)
iteration = 1;              % first iteration
mode_dominance_min = 100;   % minimum mode dominance (mode_frequency/samples_in_distribution as a percentage recommend >20)
mode_dominance = 0;         % initial mode dominance
warm_up = 20;               % number of iterations before mode is calculated (recommend >20)
%phase = zeros(iter_max:1);
%nsamples = zeros(iter_max:1);
%row_add=60;
%col_add=60;
    
phase_locator_zero=zeros(nt,nrec);
phase_locator_tuned=zeros(nt,nrec);
phase_locator_all=ones(nt,nrec);
phase_matrix=zeros(nt,nrec);

% Loop phase calculation using samples from input Segy until criteria met
tic
while ((mode_dominance<mode_dominance_min)&&(iteration<=iter_max))
    
    row_a = round(1+(nt-1)*rand);
    row_b = round(1+(nt-1)*rand);
    col_a = round(1+(nrec-1)*rand);
    col_b = round(1+(nrec-1)*rand);
    
    if row_a < row_b
        row_min = row_a;
        row_max = row_b;
    else
        row_min = row_b;
        row_max = row_a;
    end
    
    if col_a < col_b
        col_min = col_a;
        col_max = col_b;
    else
        col_min = col_b;
        col_max = col_a;
    end
    
    ndata = (row_max-row_min)*(col_max-col_min);             % sample size
    
    % Use sample if large enough
    if (row_max-row_min>nw)&&(ndata>=ndata_min)&&(ndata<=ndata_max)
        fprintf('Iteration: %d\n',iteration);                                       % display current iteration number
        sample = r_est_damped((row_min:row_max),(col_min:col_max));                 % draw sample
        sample_hilbert = r_est_damped_hilbert((row_min:row_max),(col_min:col_max)); % draw hilbert sample
        
        % Do phase rotation over all angles
        for iphase=1:nphase                                                     % loop over all phases
            phase_rad = ((iphase-1)*dphase-90)*pi/180.;                         % convert to radians and center around zero
            sample_rot = cos(phase_rad)*sample+sin(phase_rad)*sample_hilbert;   % apply phase rotation to sample
            sample_kurt(iphase)=mean(kurtosis(sample_rot));               % calculate average kurtosis at each phase
        end

        % extract max kurtosis -> gives phase
        [kurt,iphase] = max(sample_kurt);
        w0_phase = (iphase-1)*dphase-90;
        
        fprintf('Phase: %d\n',w0_phase);                                % display most recent phase estimate
        phase(iteration,1)=w0_phase;                                    % store phase estimate in vector
        nsamples(iteration,1)=ndata;                                    % store number of samples used for phase estimate
        if iteration>warm_up                                            % calculate mode if warm up completed
            [mode_phase,frequency]=mode(phase);                         % store mode and its frequency
            fprintf('Modal Phase: %d\n',mode_phase);                    % display modal phase
            fprintf('Frequency: %d\n',frequency);                       % display mode count frequency
            mode_dominance=(frequency/iteration)*100;                   % update mode dominance
            fprintf('Mode Dominance: %f\n',mode_dominance);             % display mode dominance
            fprintf('Mean Phase: %d\n',round(mean(phase)));             % display mode dominance
            secs=(toc/iteration)*(iter_max-iteration);
            hours=floor(secs/3600);
            mins=floor(((secs/3600)-floor(secs/3600))*60);
            secs=round(secs-mins*60-hours*3600);
            fprintf('Estimated time remaining (h:m:s): %d:%d:%d\n',hours,mins,secs);
        end
        
        %row_midpoint=(round((row_max-row_min)/2))+row_min;
        %col_midpoint=(round((col_max-col_min)/2))+col_min;
        %row_add_min=row_midpoint-row_add;
        %row_add_max=row_midpoint+row_add;
        %col_add_min=col_min;
        %col_add_max=col_max;
          
        phase_locator_all((row_min:row_max),(col_min:col_max))=phase_locator_all((row_min:row_max),(col_min:col_max))+1;
        %phase_locator_all((row_add_min:row_add_max),(col_add_min:col_add_max))=phase_locator_all((row_add_min:row_add_max),(col_add_min:col_add_max))+1;
        %if (w0_phase>=-15)&&(w0_phase<=15)
            %phase_locator_zero((row_min:row_max),(col_min:col_max))=phase_locator_zero((row_min:row_max),(col_min:col_max))+1;
            %phase_locator_zero((row_add_min:row_add_max),(col_add_min:col_add_max))=phase_locator_zero((row_add_min:row_add_max),(col_add_min:col_add_max))+1;
        %end
        %if (w0_phase>=-90)&&(w0_phase<=-60)
            %phase_locator_tuned((row_min:row_max),(col_min:col_max))=phase_locator_tuned((row_min:row_max),(col_min:col_max))+1;
            %phase_locator_tuned((row_add_min:row_add_max),(col_add_min:col_add_max))=phase_locator_tuned((row_add_min:row_add_max),(col_add_min:col_add_max))+1;
        %end
        phase_matrix((row_min:row_max),(col_min:col_max))=phase_matrix((row_min:row_max),(col_min:col_max))+w0_phase;
            
        iteration=iteration+1;                              % advance iteration
    end
end
toc

phase_matrix=phase_matrix./phase_locator_all;

% Set parameters for confidence bound calculations
window_increment=1;     % increment of window size used for confidence bound
count=0;                % initial number of samples within window
confidence_min=80;      % minimum confidence (quantity of data within window as a percentage of total data)
confidence=0;           % initial confidence
use_mode=0;             % use mode for statistics (=1) or use mean (=0)

if use_mode==1
    % Loop increases window size symmetrically about mode until minimum confidence level met
    while confidence<confidence_min                                     
        window_max=mode_phase+window_increment;                         % window above mode
        window_min=mode_phase-window_increment;                         % window below mode
        for i=1:length(phase)                                           % scan through stored phases
            if ((phase(i,1)>=window_min)&&(phase(i,1)<=window_max))     % check if phase falls within window
                count=count+1;                                          % count all phases within window
            end
        end
        confidence=(count/length(phase))*100;                           % calculate confidence value
        if confidence<confidence_min                                    % check if confidence level is high enough
            count=0;                                                    % reset count if confidence level no met
            window_increment=window_increment+1;                        % make window bigger for next iteration
        end
    end  
else
    % Loop increases window size symmetrically about mode until minimum confidence level met
    while confidence<confidence_min                                     
        window_max=round(mean(phase))+window_increment;                 % window above mean
        window_min=round(mean(phase))-window_increment;                 % window below mean
        for i=1:length(phase)                                           % scan through stored phases
            if ((phase(i,1)>=window_min)&&(phase(i,1)<=window_max))     % check if phase falls within window
                count=count+1;                                          % count all phases within window
            end
        end
        confidence=(count/length(phase))*100;                           % calculate confidence value
        if confidence<confidence_min                                    % check if confidence level is high enough
            count=0;                                                    % reset count if confidence level no met
            window_increment=window_increment+1;                        % make window bigger for next iteration
        end
    end  
end

% Parameters needed for plots
upper_confidence = [window_max,window_max];
lower_confidence = [window_min,window_min];
xaxis_nsampels = [min(nsamples),max(nsamples)];
if use_mode==1
    phase_plot = [mode_phase,mode_phase];
else
    phase_plot = [round(mean(phase)),round(mean(phase))];
end
xaxis_niter = [1,iter_max];

% Parameters needed for plotting wavelets and amplitude spectrum
wavelets = zeros(length(w0),3);
ftwavelets = zeros(length(w0),3);
phase_rad(1,1) = window_min*pi/180;
if use_mode==1
    phase_rad(2,1) = mode_phase*pi/180;
else
    phase_rad(2,1) = round(mean(phase))*pi/180;
end
phase_rad(3,1) = window_max*pi/180;

% Small loop to calculate wavelets with confidence bounds and ampitude spectrum
for w=1:3
    wavelets(:,w) = cos(-phase_rad(w,1))*w0+sin(-phase_rad(w,1))*imag(hilbert(w0)); % phase advance => -phase_rad
    ftwavelets(:,w) = abs(fft(wavelets(:,w)));
end

% Code to plot graphs
figure(1)
subplot(2,1,1)
plot(xaxis_nsampels,phase_plot,'LineStyle','-','Color',[1 0.4 0],'LineWidth',1.5);
hold all
plot(xaxis_nsampels,upper_confidence,'LineStyle','--','Color',[1 0.4 0]);
plot(xaxis_nsampels,lower_confidence,'LineStyle','--','Color',[1 0.4 0]);
plot(nsamples,phase,'MarkerSize',8,'MarkerFaceColor',[0.6 0.2 0],'MarkerEdgeColor',[0.6 0.2 0],'Marker','.','LineStyle','none');
title('Phase Variation with Sample Size');
xlabel('Sample Size');
ylabel('Phase (degrees)');
ylim([-100 100]);
axis tight
subplot(2,1,2)
plot(xaxis_niter,phase_plot,'LineStyle','-','Color',[0.6 0.2 0],'LineWidth',1.5);
hold all
plot(xaxis_niter,upper_confidence,'LineStyle','--','Color',[0.6 0.2 0]);
plot(xaxis_niter,lower_confidence,'LineStyle','--','Color',[0.6 0.2 0]);
plot(phase,'MarkerSize',8,'MarkerFaceColor',[1 0.4 0],'MarkerEdgeColor',[1 0.4 0],'Marker','.','LineStyle','none');
title('Estimated Phase at each Iteration');
xlabel('Iteration Number');
ylabel('Phase (degrees)');
ylim([-100 100]);
axis tight
print('-r300','-dtiffn','KPERandomDist.tiff')
%saveas(gcf,'KPERandomDist.pdf')

figure(2)
subplot(2,1,1)
plot((-nw_h+1:nw_h-1)*dt,wavelets(:,1),'LineStyle','--','Color',[0.6 0.2 0],'LineWidth',1.5)
if use_mode==1
    title('Time Domain Modal Wavelet with Phase Error Bounds (80% Confidence)');
else
    title('Time Domain Mean Wavelet with Phase Error Bounds (80% Confidence)');
end
xlabel('Time (seconds)')
ylabel('Amplitude')
hold all
plot((-nw_h+1:nw_h-1)*dt,wavelets(:,2),'LineStyle','-','Color',[1 0.4 0],'LineWidth',1.5)
plot((-nw_h+1:nw_h-1)*dt,wavelets(:,3),'LineStyle','--','Color',[0.6 0.2 0],'LineWidth',1.5)
axis tight
subplot(2,1,2)
plot((0:round(nw/2))*fsample/nw,ftwavelets(1:round(nw/2)+1,2),'LineStyle','-','Color',[1 0.4 0],'LineWidth',1.5)
if use_mode==1
    title('Amplitude Spectrum of Modal Wavelet');
else
    title('Amplitude Spectrum of Mean Wavelet');
end
xlabel('Frequency (Hz)')
ylabel('Amplitude')
axis tight
print('-r300','-dtiffn','KPERandomEstimatedWavelet.tiff')
%saveas(gcf,'KPERandomEstimatedWavelet.pdf')

%figure(3)
%imagesc(1:nrec,((1:nt)+smin)*dt,phase_locator_zero./phase_locator_all)
%title('Locations of 0 Degrees Phase Estimates (+/- 15 degrees)');
%xlabel('Trace Number')
%ylabel('Time (s)')
%colorbar
%print('-r300','-dtiffn','phase_location_low.tiff')
%saveas(gcf,'phase_location_low.pdf')

%figure(4)
%imagesc(1:nrec,((1:nt)+smin-1)*dt,phase_locator_tuned./phase_locator_all)
%title('Locations of -75 Degrees Phase Estimates (+/- 15 degrees)');
%xlabel('Trace Number')
%ylabel('Time (s)')
%colorbar
%print('-r300','-dtiffn','phase_location_high.tiff')
%saveas(gcf,'phase_location_high.pdf')

%figure(5)
%imagesc(1:nrec,((1:nt)+smin-1)*dt,phase_locator_all)
%title('Locations of All Phase Estmates');
%xlabel('Trace Number')
%ylabel('Time (s)')
%colorbar
%print('-r300','-dtiffn','phase_location_all.tiff')
%saveas(gcf,'phase_location_all.pdf')

figure(6)
imagesc(1:nrec,((1:nt)+smin-1)*dt,phase_matrix)
title('Average Phase');
xlabel('Trace Number')
ylabel('Time (s)')
colorbar
print('-r300','-dtiffn','phase_matrix.tiff')
%saveas(gcf,'phase_matrix.pdf')

figure(7)
X=min(phase):3:max(phase);
hist(phase,X);
histogram = findobj(gca,'Type','patch');
set(histogram,'FaceColor',[1 0.4 0])
title('Histogram of Phase Estimates')
xlabel('Phase Bins (degrees)')
ylabel('Count')
hold all
plot(phase_plot,0:frequency:frequency,'LineStyle','-','Color',[0.6 0.2 0],'LineWidth',1.5);
plot(upper_confidence,0:frequency:frequency,'LineStyle','--','Color',[0.6 0.2 0]);
plot(lower_confidence,0:frequency:frequency,'LineStyle','--','Color',[0.6 0.2 0]);
axis tight
print('-r300','-dtiffn','histogram_phase.tiff')
%saveas(gcf,'histogram_phase.pdf')

figure(8)
imagesc(1:nrec,((1:nt)+smin-1)*dt,phase_matrix)
hold all
wiggle(1:nrec,((1:nt)+smin-1)*dt,Segy,'wiggle',5000)
title('Average Phase with Input Segy Overlay');
xlabel('Trace Number')
ylabel('Time (s)')
colorbar
hold off
print('-r300','-dtiffn','phase_matrix_segy_overlay.tiff')
%saveas(gcf,'phase_matrix_segy_overlay.pdf')

figure(9)
imagesc(1:nrec,((1:nt)+smin-1)*dt,Segy)
colormap(flipud(gray))
title('Input Segy');
xlabel('Trace Number')
ylabel('Time (s)')
colorbar
print('-r300','-dtiffn','input_segy.tiff')
%saveas(gcf,'input_segy.pdf')

fid = fopen('KPERandom_stats.txt','wt');
fprintf(fid,'iterations: %d\n',iteration-1);
fprintf(fid,'modal phase: %d degrees\n',mode_phase);
fprintf(fid,'mode dominance: %f percent\n',mode_dominance);
fprintf(fid,'mean phase: %d degrees\n',round(mean(phase)));
fprintf(fid,'Wavelet \n',w0);
fprintf(fid,'%d \n',w0);
if use_mode==1
    fprintf(fid,'80 percent confidence (about mode): +-%d degrees',window_max-mode_phase);
else
    fprintf(fid,'80 percent confidence (about mean): +-%d degrees',window_max-round(mean(phase)));
end
fclose(fid);

save workspace

%close all
%clear -except dir

%cd ..

%end