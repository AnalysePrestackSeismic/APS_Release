close all
clear all

save_on = 0 % 0 saves no pictures, 1 saves all pictures

if save_on ==1 % makes a new directory to save the pictures in
    test_id = round(10000*rand);
    system(sprintf('mkdir ~/rnmo_test_%d',test_id));
    cd(sprintf('~/rnmo_test_%d',test_id));
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
tvw = 1; % switch for time varying wavelets (0 = off, 1 = on)

Q = 100;                                                % Q factor to be chosen                         
for i=1:length(angles)
tvw_Q(1,i) = exp(-(pi()*ricker_cf(1,i)/cos((angles(1,i)/180)*2*pi())/Q));     % the attenuation varies with frequency and incidence angle, with an assumption of t=1 and 0 degrees, t gets scaled larger for increasing offset              
end

tvw_att = tvw_Q;  % wavelet attenuation for time varying wavelets, = amplitude scaling for the final wavelet. need to assume a time/depth for the last trace to occur to do proper Q
IG_taper_length = 50;
coloured_noise = 1;                                     % switch for coloured noise (0 = off, 1 = on)
sn_ratio = [3,10,3];                                    % time varying signal to noise ratio (shallow,middle,deep)
cf_coloured_noise = [45 40 35 30 25 20];                % central frequency of angle-dependent coloured noise
for i=1:length(angles)
    tvw_att_n(1,i) = exp(-((pi()*cf_coloured_noise(1,i)/(cos((angles(1,i)/180)*2*pi())/Q))));      % Q for noise, as before. follows with the assumption that the noise is coherent, and attenuated. rather than random and uniform across angles
end

res_mov = 1;                                            % this is the switch for residual moveout, 1 is on, 0 is off

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
        tp=1/cos((angles(1,ii)/180)*2*pi());
        wp=ricker_cf(ii);
        ricker_cf_att=(wp.^2)*(sqrt((((pi()*tp)/(4*Q)).^2)+(1/wp.^2))-((pi()*tp)/(4*Q)));
        w_tmp = {ricker(ricker_cf(ii),ricker_sampling),tvw_att(ii)*ricker(ricker_cf_att,ricker_sampling)};   %has wavelet for angle set stored as original and attenuated in array cells 1 and 2
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
            tp=1/cos((angles(1,ii)/180)*2*pi());
            wp=ricker_cf(ii);
            noise_cf_att=(wp.^2)*(sqrt((((pi()*tp)/(4*Q)).^2)+(1/wp.^2))-((pi()*tp)/(4*Q)));
            w_tmp = {ricker(cf_coloured_noise(ii),ricker_sampling),tvw_att_n(ii)*ricker(noise_cf_att,ricker_sampling)}; %w_tmp is array of ricker wavelet, and attenuated ricker, not random
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
    
%Spectral Analysis of data without noise
for j=0:5            % these 2 loops take the synthetic trace from the 1 column the 6 angle gathers were in, and output into a 6 column matrix, so as to make the amplitude spectra
    for i=1:nt
syn_traces(i,j+1) = syn((i+j*nt),1);
    end
end

    syn_noise = noise_matrix*rand(2*nt,1);      %noise matrix was previously wavelets of the central frequency given in the variables, now * random numbers, different ricker wavelets for angles, proportioned above, but same random numbers to proportion them
    syn_noise = (syn_noise./repmat(sqrt(interp1([1 round(nt/2) nt],sn_ratio,(1:1:nt)','linear')),length(angles),1))*(std(syn)/std(syn_noise)); % here is noise trace which is superimposed over whole data, 6 x 1000 samples long

%Frequency Spectra of noise wavelets, after randomisation and energy
%balancing

    for j=0:5
    for i=1:nt
syn_n_traces(i,j+1) = syn_noise((i+j*nt),1);  %reads the t=0 column of each of the 6 angle sets, can be used to read out different times, puts into syn_noise. noise wavelet after randomising
    end
    end
    
syn = reshape(syn+syn_noise,nt,[]);             %noise superimposed over synthetic
else
    syn = reshape(w_matrix*[I;G],nt,[]); 
end

IG = I.*G;
IG(IG<0) = 0;
for z = 0:1
%% residual moveout simulation
if z==1
if res_mov==1
s_grid = repmat((1:1:nt)',1,6);
s_grid_i = s_grid;
s_grid_i(1:200,6) = (5:1:204)'; %moves up
s_grid_i(500:700,6) = (496:1:696)'; %moves down
s_grid_i(850:900,6) = (848:1:898);

s_grid_i(1:200,5) = (3:1:202)';
s_grid_i(500:700,5) = (498:1:698)';
s_grid_i(850:900,5) = (849:1:899)';

for ii=1:6
tmp(:,ii) = interp1(s_grid(:,ii),syn(:,ii),s_grid_i(:,ii));
end
syn = tmp;
else
end
figure (7)
imagesc(syn)

else
end

%% solve non-conditioned AVA

NCava = [ones(size(angles')) sin(angles'*pi/180).*sin(angles'*pi/180)]\syn';

NCI = NCava(1,:)';
NCG = NCava(2,:)';
NCIG = NCI.*NCG;
NCIG(NCIG<0) = 0;

%% solve conditioned AVA

ft_taper_length = 10;                                                      % assymetric (tapers only the high frequency energy to zero at Nyquist frequency)
wavelet_taper_length = 30;                                                 % symmetric about zero (outside +/- wavelet taper length about zero, amplitude is set to zero, while inside there is a cosine taper)
wavelet_length_C = 151;
half_wavelet_length = floor(wavelet_length_C/2);
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
for j=1:6   %power spectrum calculated later, this just makes the matrix of wave estimation
w_est_ava=w_est;
end

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
    wavelet_length_inv = 301;                       % variable / assumption ? or long enough to contain all viable wavelet estimations?
    half_wavelet_length = floor(wavelet_length_inv/2);
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
 
    %building matrices for inversion
    for ii = 1:length(angles)
        tvw_tmp = interp1(wavelet_grid',w_est(1+(ii-1)*wavelet_length_inv:ii*wavelet_length_inv,:)',(1:1:nt)','linear','extrap');
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
%% stuff for plots
%used in several figures, starting with (1)
scaledNCI = NCI*norm(I)/norm(NCI(subsample));
scaledNCG = NCG*norm(G)/norm(NCG(subsample));
scaledCI = CI*norm(I)/norm(CI(subsample));
scaledCG = CG*norm(G)/norm(CG(subsample));
scaledinvI = invI*norm(I)/norm(invI(subsample));
scaledinvG = invG*norm(G)/norm(invG(subsample));

% cross correlations between IG estimates and the real data, to be added
% into the titles, used in figure 2
%need to convolve I and G with appropriate wavelet estimations
Incw=I;
Gncw=G;

Icw=conv(I,w_est_ava(:,1),'same');
Gcw=conv(G,w_est_ava(:,1),'same');
Iinvw=conv(I,w_est_tr(:,1),'same');
Ginvw=conv(G,w_est_tr(:,1),'same');

corr_I_NC=max(xcorr(Incw,scaledNCI,'coeff')); % I correlation with raw amplitude regression
corr_G_NC=max(xcorr(Gncw,scaledNCG,'coeff')); % G correlation with raw regression
corr_I_C=max(xcorr(Icw,scaledCI,'coeff')); % I corr preconditioned regression
corr_G_C=max(xcorr(Gcw,scaledCG,'coeff')); % G corr preconditioned regression
corr_I_inv=max(xcorr(Iinvw,scaledinvI,'coeff')); % I corr inv
corr_G_inv=max(xcorr(Ginvw,scaledinvG,'coeff')); % G corr inv

% residual plot analysis, use original wavelets and I and G estimates to
% remake synthetic to be used for comparison
NCsyn = reshape(w_matrix(1:nt*length(angles),:)*[NCI;NCG],nt,[]);
Csyn = reshape(w_matrix(1:nt*length(angles),:)*[CI;CG],nt,[]);
invsyn = reshape(w_matrix(1:nt*length(angles),:)*[invI;invG],nt,[]);


% used in frequency spectra plot 7
n_freq = 1/(2*ricker_sampling);
f_axis = (0:n_freq/(floor(nt/2)-1):n_freq);

%power spectra of the original synthetic, the noise spectrum, the
%conditioned ava spectrum, and the inverse spectrum, all with 6 columns of
%angle
power_syn =   20*log10(bsxfun(@rdivide,abs(fft(syn_traces)),max(abs(fft(syn_traces)))));
power_syn = power_syn(1:length(f_axis),:);
power_noise = 20*log10(bsxfun(@rdivide,abs(fft(syn_n_traces)),max(abs(fft(syn_n_traces)))));
power_noise = power_noise(1:length(f_axis),:);
power_Csyn = 20*log10(bsxfun(@rdivide,abs(fft(Csyn)),max(abs(fft(Csyn)))));
power_Csyn = power_Csyn(1:length(f_axis),:);
power_invsyn = 20*log10(bsxfun(@rdivide,abs(fft(invsyn)),max(abs(fft(invsyn)))));
power_invsyn = power_invsyn(1:length(f_axis),:);


%wavelet power spectra

power_real_w = 20*log10(bsxfun(@rdivide,abs(fft(real_w)),max(abs(fft(real_w)))));
f_rw = (0:n_freq/(floor(length(real_w)/2)-1):n_freq); %axis for plots
power_w_est_C = 20*log10(bsxfun(@rdivide,abs(fft(w_est_ava)),max(abs(fft(w_est_ava)))));
f_c= (0:n_freq/((wavelet_length_C/2)-1):n_freq);
power_w_est_inv = 20*log10(bsxfun(@rdivide,abs(fft(w_est_tr)),max(abs(fft(w_est_tr)))));
f_inv= (0:n_freq/((wavelet_length_inv/2)-1):n_freq);
%% plots

axis_number_fontsize = 14;
axis_title_fontsize = 14;
title_fontsize = 16;
line_width = 2;
%% syn plots

figure(1+z)

set(1+z,'Units','inches','Position',[0 0 28 10]);

subplot(3,1,1)
plot(syn,'LineWidth',line_width)
axis tight
title('Input Synthetic Angle Traces','FontSize',title_fontsize)
xlabel('Sample number','FontSize',axis_title_fontsize)
ylabel('Amplitude','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

subplot(3,1,2)
stem(I,'marker','none','LineWidth',line_width)
axis tight
title('Input Synthetic Intercept','FontSize',title_fontsize)
xlabel('Sample number','FontSize',axis_title_fontsize)
ylabel('Amplitude','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

subplot(3,1,3)
stem(G,'marker','none','LineWidth',line_width)
axis tight
title('Input Synthetic Gradient','FontSize',title_fontsize)
xlabel('Sample number','FontSize',axis_title_fontsize)
ylabel('Amplitude','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

if save_on == 1
    set(gcf,'PaperPositionMode','auto');
    savefig('synthetic','eps','pdf','tiff','jpeg');
    hgsave('synthetic','-v7.3');
end

figure(3+z)

set(3+z,'Units','inches','Position',[0 0 23 16]);

bpI = w_matrix(1:nt,1:nt)*I;
bpG = w_matrix(1:nt,1:nt)*G;
bpI = bpI*norm(I)/norm(bpI(subsample));
bpG = bpG*norm(G)/norm(bpG(subsample));

subplot(2,2,1)
time = (1:1:nt);
% time(IG==0) = 0;
scatter(bpI(subsample),bpG(subsample),[],time(subsample),'filled','SizeData',50);
axis([-0.3 0.3 -0.3 0.3]);
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Sample Number','FontSize',axis_title_fontsize);
title('True Intercept vs Gradient filtered to seismic bandwidth','FontSize',title_fontsize)
xlabel('Intercept','FontSize',axis_title_fontsize)
ylabel('Gradient','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

subplot(2,2,2)
% time = (1:1:nt);
% time(NCIG==0) = 0;
% scatter(NCI(subsample),NCG(subsample),[],time(subsample),'filled'); axis equal
colours = ((scaledNCI - bpI).^2 + (scaledNCG - bpG).^2).^(0.5);
scatter(scaledNCI(subsample),scaledNCG(subsample),[],colours(subsample),'filled','SizeData',50);
axis([-0.3 0.3 -0.3 0.3]);
caxis([0 0.25]);
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Euclidean distance from true IG location','FontSize',axis_title_fontsize);
title({'Intercept vs Gradient from raw amplitude regression';sprintf('Number of points >%.3f from true IG location = %d',0.1,sum(colours>0.1));sprintf('Furthest point from true IG location = %.3f',max(colours))},'FontSize',title_fontsize)
xlabel('Intercept','FontSize',axis_title_fontsize)
ylabel('Gradient','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

subplot(2,2,3)
% time = (1:1:nt);
% time(CIG==0) = 0;
% scatter(CI(subsample),CG(subsample),[],time(subsample),'filled'); axis equal
colours = ((scaledCI - bpI).^2 + (scaledCG - bpG).^2).^(0.5);
scatter(scaledCI(subsample),scaledCG(subsample),[],colours(subsample),'filled','SizeData',50);
axis([-0.3 0.3 -0.3 0.3]);
caxis([0 0.25]);
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Euclidean distance from true IG location','FontSize',axis_title_fontsize);
title({'Intercept vs Gradient from conditioned amplitude regression';sprintf('Number of points >%.3f from true IG location = %d',0.1,sum(colours>0.1));sprintf('Furthest point from true IG location = %.3f',max(colours))},'FontSize',title_fontsize)
xlabel('Intercept','FontSize',axis_title_fontsize)
ylabel('Gradient','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

subplot(2,2,4)
% time = (1:1:nt);
% time(invIG==0) = 0;
% scatter(invI(subsample),invG(subsample),[],time(subsample),'filled'); axis equal
colours = ((scaledinvI - bpI).^2 + (scaledinvG - bpG).^2).^(0.5);
scatter(scaledinvI(subsample),scaledinvG(subsample),[],colours(subsample),'filled','SizeData',50);
axis([-0.3 0.3 -0.3 0.3]);
caxis([0 0.25]);
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Euclidean distance from true IG location','FontSize',axis_title_fontsize);
title({'Intercept vs Gradient from inversion';sprintf('Number of points >%.3f from true IG location = %d',0.1,sum(colours>0.1));sprintf('Furthest point from true IG location = %.3f',max(colours))},'FontSize',title_fontsize)
xlabel('Intercept','FontSize',axis_title_fontsize)
ylabel('Gradient','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

if save_on == 1
    set(gcf,'PaperPositionMode','auto')
    savefig('crossplot','eps','pdf','tiff','jpeg');
    hgsave('crossplot','-v7.3');
end

figure(5+z)

set(5+z,'Units','inches','Position',[0 0 23 16]);

ncolours = 64;
colours = [[ones(ncolours/2,1);(1:-1/((ncolours/2)-1):0)'] [(0:1/((ncolours/2)-2):1)';1;1;(1:-1/((ncolours/2)-2):0)'] [(0:1/((ncolours/2)-1):1)';ones(ncolours/2,1)]];

subplot(3,3,1)
imagesc(angles,(1:1:nt),syn,[-0.9*max(max(abs(syn))) 0.9*max(max(abs(syn)))])
axis fill
colormap(colours);
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Amplitude','FontSize',axis_title_fontsize);
title('(1) Input Synthetic Angle Traces','FontSize',title_fontsize)
xlabel('Angles (degrees)','FontSize',axis_title_fontsize)
ylabel('Sample Number','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

subplot(3,3,2)
imagesc(angles,(1:1:nt),NCsyn*norm(syn)/norm(NCsyn),[-0.9*max(max(abs(syn))) 0.9*max(max(abs(syn)))])
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Amplitude','FontSize',axis_title_fontsize);
title({'(2) Angle Traces forward modelled using';'I and G from raw amplitude regression'},'FontSize',title_fontsize)
xlabel('Angles (degrees)','FontSize',axis_title_fontsize)
ylabel('Sample Number','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

subplot(3,3,3)
imagesc(angles,(1:1:nt),syn - (NCsyn*norm(syn)/norm(NCsyn)),[-0.9*max(max(abs(syn))) 0.9*max(max(abs(syn)))])
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Amplitude','FontSize',axis_title_fontsize);
title(sprintf('Residual of (1) - (2) l_2-norm = %.2f',norm(syn - (NCsyn*norm(syn)/norm(NCsyn)))),'FontSize',title_fontsize)
xlabel('Angles (degrees)','FontSize',axis_title_fontsize)
ylabel('Sample Number','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

subplot(3,3,4)
imagesc(angles,(1:1:nt),syn,[-0.9*max(max(abs(syn))) 0.9*max(max(abs(syn)))])
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Amplitude','FontSize',axis_title_fontsize);
title('(1) Input Synthetic Angle Traces','FontSize',title_fontsize)
xlabel('Angles (degrees)','FontSize',axis_title_fontsize)
ylabel('Sample Number','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

subplot(3,3,5)
imagesc(angles,(1:1:nt),Csyn*norm(syn)/norm(Csyn),[-0.9*max(max(abs(syn))) 0.9*max(max(abs(syn)))])
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Amplitude','FontSize',axis_title_fontsize);
title({'(3) Angle Traces forward modelled using';'I and G from conditioned amplitude regression'},'FontSize',title_fontsize)
xlabel('Angles (degrees)','FontSize',axis_title_fontsize)
ylabel('Sample Number','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

subplot(3,3,6)
imagesc(angles,(1:1:nt),syn - (Csyn*norm(syn)/norm(Csyn)), [-0.9*max(max(abs(syn))) 0.9*max(max(abs(syn)))])
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Amplitude','FontSize',axis_title_fontsize);
title(sprintf('Residual of (1) - (3) l_2-norm = %.2f',norm(syn - (Csyn*norm(syn)/norm(Csyn)))),'FontSize',title_fontsize)
xlabel('Angles (degrees)','FontSize',axis_title_fontsize)
ylabel('Sample Number','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

subplot(3,3,7)
imagesc(angles,(1:1:nt),syn,[-0.9*max(max(abs(syn))) 0.9*max(max(abs(syn)))])
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Amplitude','FontSize',axis_title_fontsize);
title('(1) Input Synthetic Angle Traces','FontSize',title_fontsize)
xlabel('Angles (degrees)','FontSize',axis_title_fontsize)
ylabel('Sample Number','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

subplot(3,3,8)
imagesc(angles,(1:1:nt),invsyn*norm(syn)/norm(invsyn),[-0.9*max(max(abs(syn))) 0.9*max(max(abs(syn)))])
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Amplitude','FontSize',axis_title_fontsize);
title('(4) Angle Traces forward modelled using I and G from inversion','FontSize',title_fontsize)
xlabel('Angles (degrees)','FontSize',axis_title_fontsize)
ylabel('Sample Number','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

subplot(3,3,9)
imagesc(angles,(1:1:nt),syn - (invsyn*norm(syn)/norm(invsyn)), [-0.9*max(max(abs(syn))) 0.9*max(max(abs(syn)))])
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Amplitude','FontSize',axis_title_fontsize);
title(sprintf('Residual of (1) - (4) l_2-norm = %.2f',norm(syn - (invsyn*norm(syn)/norm(invsyn)))),'FontSize',title_fontsize)
xlabel('Angles (degrees)','FontSize',axis_title_fontsize)
ylabel('Sample Number','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

if save_on == 1
    set(gcf,'PaperPositionMode','auto')
    savefig('gathers','eps','pdf','tiff','jpeg');
    hgsave('gathers','-v7.3');

    save('workspace.mat','-v7.3');
end

end