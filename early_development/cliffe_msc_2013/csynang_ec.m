close all
clear all

save_on = 0 % 0 saves no pictures, 1 saves all pictures

if save_on ==1 % makes a new directory to save the pictures in
    test_id = round(10000*rand);
    system(sprintf('mkdir ~/digi_test_%d',test_id));
    cd(sprintf('~/digi_test_%d',test_id));
end
%% make synthetic

nt = 1000;                                              % trace length in samples
sparsity = 0.7;                                         % proportion of non-zero amplitudes (0-1) 1 = fully sparse
max_amplitude = 0.15;
ricker_cf = [30 22 20 19 17 15];                        % central frequency of ricker wavelets to use to make synthetic traces (one per angle)
ricker_sampling = 0.004;                                % time sampling rate
background_dominace = 0.8;                              % proportion of the data that sits on the background IG trend (0-1) 1 = all background
chi = 21-(0:ricker_sampling:(nt-1)*ricker_sampling)';   % background trend angle (varies with time)
chinoise = 0.8;                                           % add noise to the chi angle (0-1) 0 = no noise, spread from the trends
IGnoise = 0.8;                                           % add noise to the IG data (0-1) 0 = no noise, general spread across cross plot
angles = [10 15 20 25 30 35];                           % angle of the traces
subsample = [1;zeros((0.004/ricker_sampling)-1,1)];
subsample = bsxfun(@times,ones(1,floor(nt/length(subsample))),subsample);
subsample = logical(subsample(:));
anomaly_multiplier = 1.5;                               % boost the off-trend amplitudes (>0)
tvw = 1;                                                % switch for time varying wavelets (0 = off, 1 = on)
Q = 150;                                                % Q factor to be chosen                         
for i=1:length(angles)
tvw_Q(1,i) = exp(-(((pi()*ricker_cf(1,i)/cos(angles(1,i)/2*pi()*180)))/Q));     % the attenuation varies with frequency and incidence angle, with an assumption of t=1 and 0 degrees, t gets scaled larger for increasing offset              
end
tvw_att = tvw_Q;                                        % wavelet attenuation for time varying wavelets, = amplitude scaling for the final wavelet. need to assume a time/depth for the last trace to occur to do proper Q
IG_taper_length = 50;
coloured_noise = 1;                                     % switch for coloured noise (0 = off, 1 = on)
sn_ratio = [3,10,3];                                    % time varying signal to noise ratio (shallow,middle,deep)
cf_coloured_noise = [45 40 35 30 25 20];                % central frequency of angle-dependent coloured noise
for i=1:length(angles)
    tvw_att_n(1,i) = exp(-((pi()*cf_coloured_noise(1,i)/cos(angles(1,i)/2*pi()*180))/Q));      % Q for noise, as before. follows with the assumption that the noise is coherent, and attenuated. rather than random and uniform across angles
end
ref = synref(nt,sparsity,max_amplitude);                % make reflection coefficients
class = nan(nt,1);                                      % class starts as a series of nt nans
random = rand(nt,1);                                    % series of nt random numbers
class(ref~=0) = random(ref~=0);                         % where ref isn't zero, insert corresponding random number
class(class>=(1-background_dominace)) = 1;              % replace values in class vector with 1 or 2 depending on background dominance, ensures the right proportion of anomalous points
class(class<(1-background_dominace)) = 2;               
class_2_flag = 0;
for ii = 1:nt                                           %this loop makes anomalies, those with class 2 are flagged wtih a number 1, reflectivity is made negative, anomalous and squared
    if class(ii) == 2
        class_2_flag = 1;
        ref(ii) = -anomaly_multiplier*abs(ref(ii));
    end
    if class(ii) == 1                                   % class 1, 
        if class_2_flag == 1                            % if the previous loop has just made an anomalous value, the next one has class made =3, +ve anomalous and squared
            class(ii) = 3;
            ref(ii) = anomaly_multiplier*abs(ref(ii));  % ensures equal number of 2s and 3s, and equal +ve and -ve - white earth reflectivity assumption, could bias this at this point, break assumption
            class_2_flag = 0;                           % flag reset to zero and looped
        end
    end
end

I = ref;                                                %intercept is zero offset reflectivity, ie ref
G = zeros(nt,1);  
chi_array = chi+(chinoise-2*chinoise*rand(nt,1));                    %chi angles for each point in trace length, chi angle + noise, +random factor scaled by noise value, different for every point
G(class==1) = -I(class==1)./tan(chi_array(class==1)*pi/180);         %converts chi angle for each point to a gradient depending on the class. class 1 is background, 2 is negative anomaly, 3 is positive anomaly, 
G(class==2) = -I(class==2)./tan(90+chi_array(class==2)*pi/180);      % anomaly is perpendicular to background
G(class==3) = -I(class==3)./tan(90+chi_array(class==3)*pi/180);      %G therefore scattered by chinoise, not equal to chi for background

Gnoise = IGnoise*std(G)*randn(nt,1);                                 % gradient noise vector, std proportions, as does IGnoise value, and random extra proportion
G(~isnan(class)) = G(~isnan(class))+Gnoise(~isnan(class));           %adds noise onto the G vector at all the points where reflectivity occurs

Inoise = IGnoise*std(I)*randn(nt,1);
I(~isnan(class)) = I(~isnan(class))+Inoise(~isnan(class));           %adds noise onto the I vector at all the points where reflectivity occurs


I = bsxfun(@times,[fliplr(0.5+0.5*cos(0:pi/IG_taper_length:pi)),ones(1,nt-2*IG_taper_length-2),0.5+0.5*cos(0:pi/IG_taper_length:pi)]',I);
G = bsxfun(@times,[fliplr(0.5+0.5*cos(0:pi/IG_taper_length:pi)),ones(1,nt-2*IG_taper_length-2),0.5+0.5*cos(0:pi/IG_taper_length:pi)]',G);

eref = bsxfun(@plus,I,bsxfun(@times,G,sin(angles*pi/180).*sin(angles*pi/180)));

if tvw == 1
    for ii = 1:length(angles) %for every angle gather
        w_tmp = {ricker(ricker_cf(ii),ricker_sampling),tvw_att(ii)*ricker(tvw_att(1,ii)*ricker_cf(ii),ricker_sampling)};   %has wavelet for angle set stored as original and attenuated in array cells 1 and 2
        w_tmp{1} = circshift([w_tmp{1};zeros(length(w_tmp{2})-length(w_tmp{1}),1)],(length(w_tmp{2})-length(w_tmp{1}))/2); %puts zeros at begining and end
         ar_tmp(:,ii)=w_tmp(1,1);   % this bit saves the original wavelet which is used to make the synthetic before time varying, so it can be compared to the wavelet estimated by the methods
        ar_tmp2(:,ii)=w_tmp(1,2);
        for j=1:size(ar_tmp{1,ii})
        real_w(j,ii)=ar_tmp{1,ii}(j,1);
        end
        for j=1:size(ar_tmp2{1,ii})
        real_w_att(j,ii)=ar_tmp2{1,ii}(j,1); %saves the attenuated original wavelet too
        end
        
        w_tmp = interp1([1;nt],[w_tmp{1}';w_tmp{2}'],(1:1:nt)');                               % interpolates, between original wavelet and attenuated wavelet, fills huge matrix
        if ii == 1
            w_matrix = [spdiags(w_tmp,-floor(size(w_tmp,2)/2):floor(size(w_tmp,2)/2),nt,nt), spdiags(bsxfun(@times,sind(angles(ii))*sind(angles(ii)),w_tmp),-floor(size(w_tmp,2)/2):floor(size(w_tmp,2)/2),nt,nt)]; %manipulates w_tmp for angle set into the wavelet matrix
        else
            w_matrix = [w_matrix; [spdiags(w_tmp,-floor(size(w_tmp,2)/2):floor(size(w_tmp,2)/2),nt,nt), spdiags(bsxfun(@times,sind(angles(ii))*sind(angles(ii)),w_tmp),-floor(size(w_tmp,2)/2):floor(size(w_tmp,2)/2),nt,nt)]];
        end
    end
else
    for ii = 1:length(angles)
        w_tmp = ricker(ricker_cf(ii),ricker_sampling);
        if ii == 1
            w_matrix = [spdiags(repmat(w_tmp',nt,1),-floor(length(w_tmp)/2):floor(length(w_tmp)/2),nt,nt), spdiags(repmat(sind(angles(ii))*sind(angles(ii))*w_tmp',nt,1),-floor(length(w_tmp)/2):floor(length(w_tmp)/2),nt,nt)];
        else
            w_matrix = [w_matrix; [spdiags(repmat(w_tmp',nt,1),-floor(length(w_tmp)/2):floor(length(w_tmp)/2),nt,nt), spdiags(repmat(sind(angles(ii))*sind(angles(ii))*w_tmp',nt,1),-floor(length(w_tmp)/2):floor(length(w_tmp)/2),nt,nt)]];
        end
    end
end

if coloured_noise == 1
    if tvw == 1
        for ii = 1:length(angles)
            w_tmp = {ricker(cf_coloured_noise(ii),ricker_sampling),tvw_att_n(ii)*ricker(tvw_att_n(ii)*cf_coloured_noise(ii),ricker_sampling)}; %w_tmp is array of ricker wavelet, and attenuated ricker, not random
            w_tmp{1} = circshift([w_tmp{1};zeros(length(w_tmp{2})-length(w_tmp{1}),1)],(length(w_tmp{2})-length(w_tmp{1}))/2); %adds zeroes at each end, probably other changes
            ar_tmp(:,ii)=w_tmp(1,1);   % this bit saves the original wavelet which is used to make the synthetic before time varying, so it can be compared to the wavelet estimated by the methods
        ar_tmp2(:,ii)=w_tmp(1,2);
        for j=1:size(ar_tmp{1,ii})
        noise_w(j,ii)=ar_tmp{1,ii}(j,1);
        end
        for j=1:size(ar_tmp2{1,ii})
        noise_w_att(j,ii)=ar_tmp2{1,ii}(j,1); %saves the attenuated original wavelet too
        end
            
            
            w_tmp = interp1([1;nt],[w_tmp{1}';w_tmp{2}'],(1:1:nt)'); %interpolates, then addds into huge noise matrix for each point in next lines
            if ii == 1
                noise_matrix = [spdiags(w_tmp,-floor(size(w_tmp,2)/2):floor(size(w_tmp,2)/2),nt,nt), spdiags(bsxfun(@times,sind(angles(ii))*sind(angles(ii)),w_tmp),-floor(size(w_tmp,2)/2):floor(size(w_tmp,2)/2),nt,nt)];
            else
                noise_matrix = [noise_matrix; [spdiags(w_tmp,-floor(size(w_tmp,2)/2):floor(size(w_tmp,2)/2),nt,nt), spdiags(bsxfun(@times,sind(angles(ii))*sind(angles(ii)),w_tmp),-floor(size(w_tmp,2)/2):floor(size(w_tmp,2)/2),nt,nt)]];
            end
        end
    else
        for ii = 1:length(angles)
            w_tmp = ricker(cf_coloured_noise(ii),ricker_sampling);
            if ii == 1
                noise_matrix = [spdiags(repmat(w_tmp',nt,1),-floor(length(w_tmp)/2):floor(length(w_tmp)/2),nt,nt), spdiags(repmat(sind(angles(ii))*sind(angles(ii))*w_tmp',nt,1),-floor(length(w_tmp)/2):floor(length(w_tmp)/2),nt,nt)];
            else
                noise_matrix = [noise_matrix; [spdiags(repmat(w_tmp',nt,1),-floor(length(w_tmp)/2):floor(length(w_tmp)/2),nt,nt), spdiags(repmat(sind(angles(ii))*sind(angles(ii))*w_tmp',nt,1),-floor(length(w_tmp)/2):floor(length(w_tmp)/2),nt,nt)]];
            end
        end
    end
    syn = w_matrix*[I;G];     %synthetic trace, wavelet*intercept and gradient
    
    
%     figure;  %plotting the 6 synthetic traces with no noise, in 1 long trace
%     plot(syn)
% axis tight
% title('synthetic before noise')

%Spectral Analysis of data without noise
for j=0:5            % these 2 loops take the synthetic trace from the 1 column the 6 angle gathers were in, and output into a 6 column matrix, so as to make the amplitude spectra
    for i=1:1000
syn_traces(i,j+1) = syn((i+j*1000),1);
    end
end

for j=1:6               %here, power spectrum of each angle gather is computed plotting the 6 synthetic traces with no noise, in 1 long trace
Y=fft(syn_traces(:,j));
FS_syn(:,j)=Y.*conj(Y)/1000;
end

% for j = 1:6
% subplot(3,2,j)
% plot(1:500,PSpec_Y(1:500,j))
% end
% title('frequency spectra for the 6 angle gathers')

% Difference between frequency spectra
% figure;
% plot(1:500,PSpec_Y(1:500,1)-PSpec_Y(1:500,6))      % displaying the difference between the smallest and largest angle spectra
% title('difference between spectra of 1 and 6')

    syn_noise = noise_matrix*rand(2*nt,1);      %noise matrix was previously wavelets of the central frequency given in the variables, now * random numbers, different ricker wavelets for angles, proportioned above, but same random numbers to proportion them
    syn_noise = (syn_noise./repmat(sqrt(interp1([1 round(nt/2) nt],sn_ratio,(1:1:nt)','linear')),length(angles),1))*(std(syn)/std(syn_noise)); % here is noise trace which is superimposed over whole data, 6 x 1000 samples long

%Frequency Spectra of noise wavelets, after randomisation and energy
%balancing

    for j=0:5
    for i=1:1000
syn_n_traces(i,j+1) = syn_noise((i+j*1000),1);  %reads the t=0 column of each of the 6 angle sets, can be used to read out different times, puts into syn_noise. noise wavelet after randomising
    end
    end
    
Y=0;
    for j=1:6               %here, power spectrum of each angle gather is computed
Y=fft(syn_n_traces(:,j));
FS_noise(:,j)=Y.*conj(Y)/1000;
    end
%     for j=1:6
% subplot(3,2,j)
% plot(1:500,PSpec_Y_n(1:500,j))
%     end
% title('noise power spectrum, for each angle set')

syn = reshape(syn+syn_noise,nt,[]);             %noise superimposed over synthetic
else
    syn = reshape(w_matrix*[I;G],nt,[]); 
end

IG = I.*G;
IG(IG<0) = 0;

%% solve non-conditioned AVA

NCava = [ones(size(angles')) sin(angles'*pi/180).*sin(angles'*pi/180)]\syn';

NCI = NCava(1,:)';
NCG = NCava(2,:)';
NCIG = NCI.*NCG;
NCIG(NCIG<0) = 0;

%% solve conditioned AVA

ft_taper_length = 10;                                                      % assymetric (tapers only the high frequency energy to zero at Nyquist frequency)
wavelet_taper_length = 30;                                                 % symmetric about zero (outside +/- wavelet taper length about zero, amplitude is set to zero, while inside there is a cosine taper)
wavelet_length = 151;
half_wavelet_length = floor(wavelet_length/2);
hnt = floor(nt/2);
ft_taper = [ones(hnt-ft_taper_length,1); ((1+cos((0:1/(ft_taper_length-1):1)*pi)')/2)];
ft_taper = [ft_taper; zeros(nt-length(ft_taper),1)];

ft_traces_tmp = abs(fft(syn));                                             %syn is the synthetic data produced earlier, now being analysed
ft_traces_tmp = bsxfun(@times,ft_traces_tmp,ft_taper);                     % ensures zero energy at Nyquist frequency
w_est = circshift(ifft(ft_traces_tmp,'symmetric'),hnt); 
[~, peak_index] = max(w_est);
peak_index = max(peak_index);
w_est = w_est(peak_index-half_wavelet_length:peak_index+half_wavelet_length,:);
wavelet_taper = [((1+cos((-1:1/(wavelet_taper_length-1):0)*pi)')/2); ones(13,1); ((1+cos((0:1/(wavelet_taper_length-1):1)*pi)')/2)];
wavelet_taper = [zeros((size(w_est,1)-length(wavelet_taper))/2,1); wavelet_taper; zeros((size(w_est,1)-length(wavelet_taper))/2,1)];
w_est = bsxfun(@times,w_est,wavelet_taper);                                % ensures zero amplitude at terminations
w_est(isnan(w_est)) = 0;
%no wavelet normalisation here, but still reduces amlplitude
for j=1:6
w_est_ava=w_est;
    Y=fft(w_est(:,j));
FS_w_est_ava(:,j)=abs(Y(1:75,1)/151);
end
%plot(1:75,FS_w_est_ava(1:75,1))

for ii = 1:length(angles)
    filter(:,ii) = match_call(w_est(:,round((length(angles))/2)),w_est(:,ii),-6); % match_call matches an input wavelet to a desired output wavelet, rms stays the same
end

for ii = 1:length(angles)
    Csyn(:,ii) = conv(syn(:,ii),filter(:,ii),'same');                       %convolve synthetic with matching filter, edits data
    fscalar = norm(syn(:,ii))/norm(Csyn(:,ii));                             %norm (vector length) of synthetic / norm of spectral matched synthetic
    filter(:,ii) = filter(:,ii)*fscalar;                                    %adjust filter by multiplying with above
    Csyn(:,ii) = conv(syn(:,ii),filter(:,ii),'same');                       % replace adjusted data from step 1 with adjusted filtered version
end

% Csyn equivalent to synthetic, 6 sets of angle gathers, already in
% columns, prior to inversion for I and G
for j=1:6
Y=fft(Csyn(:,j));
FS_Csyn(:,j)=Y.*conj(Y)/1000;
end

Cava = [ones(size(angles')) sin(angles'*pi/180).*sin(angles'*pi/180)]\Csyn'; % inversion

CI = Cava(1,:)';
CG = Cava(2,:)';
CIG = CI.*CG;
CIG(CIG<0) = 0;

%% Inverted AVA
wsmooth = 0;                                    % variable, controls smoothness, but reduces how good the fit is, never do this on real data prior to checking the unsmoothed data
n_wavelets = 3;

smooth = spdiags([-wsmooth*ones(2*nt,1) 2*wsmooth*ones(2*nt,1) -wsmooth*ones(2*nt,1)],[-1 0 1],2*nt,2*nt);

if tvw == 1                                    % switch for time varying wavelets (0 = off, 1 = on)
    wavelet_step = floor(nt/(n_wavelets+1));
    wavelet_grid = (wavelet_step:wavelet_step:n_wavelets*wavelet_step);
    ft_taper_length = 10;                       % assymetric (tapers only the high frequency energy to zero at Nyquist frequency)
    wavelet_taper_length = 30;                  % symmetric about zero (outside +/- wavelet taper length about zero, amplitude is set to zero, while inside there is a cosine taper)
    wavelet_length = 301;                       % variable / assumption ? or long enough to contain all viable wavelet estimations?
    half_wavelet_length = floor(wavelet_length/2);
    hnt = wavelet_step;
    ft_taper = [ones(hnt-ft_taper_length,1); ((1+cos((0:1/(ft_taper_length-1):1)*pi)')/2)];
    ft_taper = [ft_taper; zeros(2*wavelet_step-length(ft_taper),1)];
    
    for jj = 1:n_wavelets
        ft_traces_tmp = abs(fft(syn(1+(jj-1)*wavelet_step:(jj+1)*wavelet_step,:)));
        ft_traces_tmp = bsxfun(@times,ft_traces_tmp,ft_taper);              % ensures zero energy at Nyquist frequency
        w_est = circshift(ifft(ft_traces_tmp,'symmetric'),hnt);             % wavelet estimations, inverse FT of the abs of FT of synthetic*wavelet_step
        [~, peak_index] = max(w_est);                                       % finds maxima, output into peak index
        peak_index = max(peak_index);
        w_est = w_est(peak_index-half_wavelet_length:peak_index+half_wavelet_length,:);
        wavelet_taper = [((1+cos((-1:1/(wavelet_taper_length-1):0)*pi)')/2); ones(13,1); ((1+cos((0:1/(wavelet_taper_length-1):1)*pi)')/2)]; 
        wavelet_taper = [zeros((size(w_est,1)-length(wavelet_taper))/2,1); wavelet_taper; zeros((size(w_est,1)-length(wavelet_taper))/2,1)];   %designs taper
        w_est = bsxfun(@times,w_est,wavelet_taper);      %301x6 at this point                  % ensures zero amplitude at terminations, w_est*taper
        w_est(isnan(w_est)) = 0; 
        normalise = sum(sqrt(w_est.^2));                                   
        w_est = normalise(1)*bsxfun(@rdivide,w_est,normalise); %makes the estimated wavelet same energy as the original? leaves amplitude spectra a lot lower amplitude though
        w_est = w_est(:);     %1806x1 , 6 301 length wavelets in a line
        w_est_tmp(:,jj) = w_est;
    end
     
    w_est = w_est_tmp;
    
    
    % Calculating the difference between the original and the estimated
    % frequency spectrum
    for j=0:5
        for i=1:301
    w_est_tr(i,j+1)=w_est(i+j*301,1);
        end
    end
    for j=1:6
    Y=fft(w_est_tr(:,j)/301);
    FS_w_est_digi(:,j)=abs(Y);
    end
    
    %plot(1:151,FS_w_est_digi(1:151,1))
    
    %building matrices for inversion
    for ii = 1:length(angles)
        tvw_tmp = interp1(wavelet_grid',w_est(1+(ii-1)*wavelet_length:ii*wavelet_length,:)',(1:1:nt)','linear','extrap');
        if ii == 1
            IGmatrix = [spdiags(tvw_tmp,-floor(size(tvw_tmp,2)/2):floor(size(tvw_tmp,2)/2),nt,nt), spdiags(bsxfun(@times,sind(angles(ii))*sind(angles(ii)),tvw_tmp),-floor(size(tvw_tmp,2)/2):floor(size(tvw_tmp,2)/2),nt,nt)];
        else
            IGmatrix = [IGmatrix; [spdiags(tvw_tmp,-floor(size(tvw_tmp,2)/2):floor(size(tvw_tmp,2)/2),nt,nt), spdiags(bsxfun(@times,sind(angles(ii))*sind(angles(ii)),tvw_tmp),-floor(size(tvw_tmp,2)/2):floor(size(tvw_tmp,2)/2),nt,nt)]];
        end
    end
    IGmatrix = [IGmatrix; smooth];
else
    for ii = 1:length(angles)
        w_tmp = w_est(:,ii);
        if ii == 1
            IGmatrix = [spdiags(repmat(w_tmp',nt,1),-floor(length(w_tmp)/2):floor(length(w_tmp)/2),nt,nt), spdiags(repmat(sind(angles(ii))*sind(angles(ii))*w_tmp',nt,1),-floor(length(w_tmp)/2):floor(length(w_tmp)/2),nt,nt)];
        else
            IGmatrix = [IGmatrix; [spdiags(repmat(w_tmp',nt,1),-floor(length(w_tmp)/2):floor(length(w_tmp)/2),nt,nt), spdiags(repmat(sind(angles(ii))*sind(angles(ii))*w_tmp',nt,1),-floor(length(w_tmp)/2):floor(length(w_tmp)/2),nt,nt)]];
        end
    end
    IGmatrix = [IGmatrix; smooth];
end

data = [syn(:);zeros(2*nt,1)];
model = [NCI;NCG];

invava = lsqr(IGmatrix,data,1e-3,500,[],[],model);  %AVA inversion, uses unconditioned AVA inversion as starting model, and synthetic as input data

invI = invava(1:nt); 
invG = invava(nt+1:end);

invIG = invI.*invG;
invIG(invIG<0) = 0;

