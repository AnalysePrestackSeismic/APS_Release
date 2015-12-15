% make sure rand_param has been run once, and use the saved variables for
% future tests, holding the ref and rand constant
clear all
close all
load param.mat

save_on = 1 % 0 saves no pictures, 1 saves all pictures

if save_on ==1 % makes a new directory to save the pictures in
    cd('~/key plots');
end
%% make synthetic
param={[30,30,30],[20,20,20],[10,10,10],[5,5,5],[2.5,2.5,2.5],[2.25,2.25,2.25],[0.75,0.75,0.75],[0.5,0.5,0.5],[0.4,0.4,0.4],[0.2,0.2,0.2],[0.1,0.1,0.1],[0.07,0.07,0.07]};
for tvw=0:1
for zz=1:length(param)

ricker_cf = [30 22 20 19 17 15];                        
ricker_sampling = 0.004;                                                    
chi = 21-(0:ricker_sampling:(nt-1)*ricker_sampling)';
chinoise = 0;
IGnoise = 0; 
w_noise=0;
w_noise_amp=0;
angles = [10 15 20 25 30 35];
Q = 150;      
cf_coloured_noise = [45 40 35 30 25 20];  
coloured_noise = 1; 
IG_taper_length = 50;
anomaly_multiplier = 1.5; 
res_mov = 0;
multiples=0;
sn_ratio = param{1,zz};
pre_white=1;
M_a=1;
ph=[0,0,0,0,0,0];                                  %phase shift in degrees for each subsequent angle, phase shift dependent on angle, doesn't vary with depth
subsample = [1;zeros((0.004/ricker_sampling)-1,1)];
subsample = bsxfun(@times,ones(1,floor(nt/length(subsample))),subsample);
subsample = logical(subsample(:));

for i=1:length(angles)
tvw_Q(1,i) = exp(-(pi()*ricker_cf(1,i)/cos((angles(1,i)/180)*2*pi())/Q));     % the attenuation varies with frequency and incidence angle, with an assumption of t=1 and 0 degrees, t gets scaled larger for increasing offset              
end

tvw_att = tvw_Q;  % wavelet attenuation for time varying wavelets, = amplitude scaling for the final wavelet. need to assume a time/depth for the last trace to occur to do proper Q
                                    % switch for coloured noise (0 = off, 1 = on)

                  % central frequency of angle-dependent coloured noise
for i=1:length(angles)
    tvw_att_n(1,i) = exp(-((pi()*cf_coloured_noise(1,i)/(cos((angles(1,i)/180)*2*pi())/Q))));      % Q for noise, as before. follows with the assumption that the noise is coherent, and attenuated. rather than random and uniform across angles
end


         


I = ref;                                                %intercept is zero offset reflectivity, ie ref
G = zeros(nt,1);  
chi_array = chi+(chinoise-2*chinoise*R1);                    %chi angles for each point in trace length, chi angle + noise, +random factor scaled by noise value, different for every point
G(class==1) = -I(class==1)./tan(chi_array(class==1)*pi/180);         %converts chi angle for each point to a gradient depending on the class. class 1 is background, 2 is negative anomaly, 3 is positive anomaly, 
G(class==2) = -I(class==2)./tan(90+chi_array(class==2)*pi/180);      % anomaly is perpendicular to background
G(class==3) = -I(class==3)./tan(90+chi_array(class==3)*pi/180);      %G therefore scattered by chinoise, not equal to chi for background

Gnoise = IGnoise*std(G)*R2;                                 % gradient noise vector, std proportions, as does IGnoise value, and random extra proportion
G(~isnan(class)) = G(~isnan(class))+Gnoise(~isnan(class));           %adds noise onto the G vector at all the points where reflectivity occurs

Inoise = IGnoise*std(I)*R3;
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
        w_tmp{1}=w_tmp{1}*cosd(ph(ii))+imag(hilbert(w_tmp{1}))*sind(ph(ii));
        w_tmp{2}=w_tmp{2}*cosd(ph(ii))+imag(hilbert(w_tmp{2}))*sind(ph(ii));
      
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
    real_w=zeros(37,6);
    for ii = 1:length(angles)
        w_tmp = ricker(ricker_cf(ii),ricker_sampling);
        if ii == 1
            real_w(1:length(w_tmp),ii)=w_tmp;
            w_matrix = [spdiags(repmat(w_tmp',nt,1),-floor(length(w_tmp)/2):floor(length(w_tmp)/2),nt,nt), spdiags(repmat(sind(angles(ii))*sind(angles(ii))*w_tmp',nt,1),-floor(length(w_tmp)/2):floor(length(w_tmp)/2),nt,nt)];
        else
            real_w(1:length(w_tmp),ii) =w_tmp;
            w_matrix = [w_matrix; [spdiags(repmat(w_tmp',nt,1),-floor(length(w_tmp)/2):floor(length(w_tmp)/2),nt,nt), spdiags(repmat(sind(angles(ii))*sind(angles(ii))*w_tmp',nt,1),-floor(length(w_tmp)/2):floor(length(w_tmp)/2),nt,nt)]];
        end
    end
end

if coloured_noise == 1
    if tvw == 1
        for ii = 1:length(angles)
            tp=1/cos((angles(1,ii)/180)*2*pi());
            wp=cf_coloured_noise(ii);
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
    %% make multiple   
    if multiples==1;
    if zz==1
        m_amp=zeros(1000,6);
    end
    
    m_amp(100+(zz-1)*100,1)=-0.1*M_a;
    m_amp(104+(zz-1)*100,2)=-0.1*M_a;
    m_amp(110+(zz-1)*100,3)=-0.12*M_a;
    m_amp(120+(zz-1)*100,4)=-0.14*M_a;
    m_amp(131+(zz-1)*100,5)=-0.16*M_a;
    m_amp(142+(zz-1)*100,6)=-0.18*M_a;
    
    tmp=zeros(1000,6);
    if tvw==1
    for i=1:6
    tmp(:,i)=conv(m_amp(:,i),real_w(:,i),'same');
    end
    mult=zeros(6000,1);
    for i=1:6
    mult((1+(i-1)*1000):1000*i,1)=tmp(:,i);
    end
    
    else
         for i=1:6
    tmp(:,i)=conv(m_amp(:,i),real_w(:,1),'same');
    end
    mult=zeros(6000,1);
    for i=1:6
    mult((1+(i-1)*1000):1000*i,1)=tmp(:,i)';
    end
    end
        
    syn = bsxfun(@plus,(w_matrix*[I;G]),mult);     %synthetic trace, wavelet*intercept and gradient
    else
       syn = w_matrix*[I;G]; 
    end
    
    
%Spectral Analysis of data without noise

    syn_noise = noise_matrix*R4;      %noise matrix was previously wavelets of the central frequency given in the variables, now * random numbers, different ricker wavelets for angles, proportioned above, but same random numbers to proportion them
    syn_noise = (syn_noise./repmat(sqrt(interp1([1 round(nt/2) nt],sn_ratio,(1:1:nt)','linear')),length(angles),1))*(std(syn)/std(syn_noise)); % here is noise trace which is superimposed over whole data, 6 x 1000 samples long

%Frequency Spectra of noise wavelets, after randomisation and energy
%balancing

    for j=0:5
    for i=1:nt
syn_n_traces(i,j+1) = syn_noise((i+j*nt),1);  %reads the t=0 column of each of the 6 angle sets, can be used to read out different times, puts into syn_noise. noise wavelet after randomising
    end
    end
    
    if w_noise==1;               %coloured noise and white noise, or just coloured noise
    syn_noise=bsxfun(@plus,syn_noise,(R5*w_noise_amp));
    syn = reshape(syn+syn_noise,nt,[]);
    
    else
    syn = reshape(syn+syn_noise,nt,[]);             %noise superimposed over synthetic
    end
        
else
    if w_noise==1;  %just white noise, or no noise
    syn_noise=R5*w_noise_amp;
    end
    if multiples==1;
       
    % first multiple
    if zz==1
        m_amp=zeros(1000,6);
    end
    m_amp(100+(zz-1)*100,1)=-0.1*M_a;
    m_amp(104+(zz-1)*100,2)=-0.1*M_a;
    m_amp(110+(zz-1)*100,3)=-0.12*M_a;
    m_amp(120+(zz-1)*100,4)=-0.14*M_a;
    m_amp(131+(zz-1)*100,5)=-0.16*M_a;
    m_amp(142+(zz-1)*100,6)=-0.18*M_a;
    tmp=zeros(1000,6);
    
 if tvw==1
    for i=1:6
    tmp(:,i)=conv(m_amp(:,i),real_w(:,i),'same');
    end
    mult=zeros(6000,1);
    for i=1:6
    mult((1+(i-1)*1000):1000*i,1)=tmp(:,i);
    end
    else
     for i=1:6
    tmp(:,i)=conv(m_amp(:,i),real_w(:,1),'same');
    end
    mult=zeros(6000,1);
    for i=1:6
    mult((1+(i-1)*1000):1000*i,1)=tmp(:,i);
    end
 end
 
 
    
if w_noise==0
    syn_noise =zeros(6000,1);
end

syn_noise = bsxfun(@plus,syn_noise,mult);   %synthetic trace, wavelet*intercept and gradient
    syn = reshape(w_matrix*[I;G]+syn_noise,nt,[]); 
    syn_traces=syn; %again, just for the amplitude spectra, must be a better way of doing this rather than dropping in this line everywhere ¬¬
  
   
    else
  
if w_noise==0;
    syn_noise =zeros(6000,1);
end
    syn = reshape(w_matrix*[I;G]+syn_noise,nt,[]); 
    syn_traces=syn; %again, just for the amplitude spectra, must be a better way of doing this rather than dropping in this line everywhere ¬¬

    end
end

IG = I.*G;
IG(IG<0) = 0;
%% residual moveout simulation
if res_mov==1
s_grid = repmat((1:1:nt)',1,6);
s_grid_i = s_grid;
s_grid_i(1:200,6) = (4:1:203)'; %moves up
s_grid_i(850:900,6) = (848:1:898);

s_grid_i(1:200,5) = (2:1:201)';
% s_grid_i(500:700,5) = (498:1:698)';
 s_grid_i(850:900,5) = (849:1:899)';

for ii=1:6
tmp(:,ii) = interp1(s_grid(:,ii),syn(:,ii),s_grid_i(:,ii));
end
syn = tmp;

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
    for j=0:length(angles)-1
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
else                                           % if time varying wavelets is turned off, no wavelet estimation is done by inv, uses wavelet estimation done by conditioned
    for ii = 1:length(angles)
        w_tmp = w_est(:,ii);
        wavelet_length_inv=wavelet_length_C;
    for j=1:length(angles)
    w_est_tr(:,j)=w_est(:,j);
    end
        if ii == 1
            IGmatrix = [spdiags(repmat(w_tmp',nt,1),-floor(length(w_tmp)/2):floor(length(w_tmp)/2),nt,nt), spdiags(repmat(sind(angles(ii))*sind(angles(ii))*w_tmp',nt,1),-floor(length(w_tmp)/2):floor(length(w_tmp)/2),nt,nt)];
        else
            IGmatrix = [IGmatrix; [spdiags(repmat(w_tmp',nt,1),-floor(length(w_tmp)/2):floor(length(w_tmp)/2),nt,nt), spdiags(repmat(sind(angles(ii))*sind(angles(ii))*w_tmp',nt,1),-floor(length(w_tmp)/2):floor(length(w_tmp)/2),nt,nt)]];
        end
    end
    IGmatrix = [IGmatrix; smooth];
end

if pre_white==1;
pw = spdiags(ones(nt,2),[0,nt],nt,2*nt);
tmp = sind(angles).*sind(angles);
tmp = repmat(tmp,nt,1);
tmp = tmp(:);
pw = repmat(pw,length(angles),1);
pw(:,1001:end) = bsxfun(@times,pw(:,1001:end),tmp);
pw = [pw;zeros(2*nt,2*nt)];
    
    IGmatrix = IGmatrix+pw;
else
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

% residual plot analysis, use original wavelets and I and G estimates to
% remake synthetic to be used for comparison
NCsyn = reshape(w_matrix(1:nt*length(angles),:)*[NCI;NCG],nt,[]);
Csyn = reshape(w_matrix(1:nt*length(angles),:)*[CI;CG],nt,[]);
invsyn = reshape(w_matrix(1:nt*length(angles),:)*[invI;invG],nt,[]);

% cross correlations between IG estimates and the real data, to be added
% into the titles, used in figure 2
%need to convolve I and G with real wavelet, to isolate the I and G
%estimating capabilities. wavelet estimating properties less important

%%convolving the original I and G values with estimated wavelets by the 2
%%methods, then cross correlation compares these to the estimated I and G.
%%should use true wavelet, so tests just the I and G estimating
%%capabilities

Incw=conv(I,real_w(:,1),'same');      %only uses first angle set of wavelets, and isn't time varying
Gncw=conv(I,real_w(:,1),'same');
Icw=conv(I,real_w(:,1),'same');
Gcw=conv(G,real_w(:,1),'same');
Iinvw=conv(I,real_w(:,1),'same');
Ginvw=conv(G,real_w(:,1),'same');

corr_I_NC=max(xcorr(Incw,scaledNCI,'coeff')); % I correlation with raw amplitude regression
corr_G_NC=max(xcorr(Gncw,scaledNCG,'coeff')); % G correlation with raw regression
corr_I_C=max(xcorr(Icw,scaledCI,'coeff')); % I corr preconditioned regression
corr_G_C=max(xcorr(Gcw,scaledCG,'coeff')); % G corr preconditioned regression
corr_I_inv=max(xcorr(Iinvw,scaledinvI,'coeff')); % I corr inv
corr_G_inv=max(xcorr(Ginvw,scaledinvG,'coeff')); % G corr inv

corr_syn_NCnr=max(xcorr(NCsyn(:,1),syn(:,1),'coeff')); % I correlation with raw amplitude regression
corr_syn_Cnr=max(xcorr(Csyn(:,1),syn(:,1),'coeff')); % I corr preconditioned regression
corr_syn_invnr=max(xcorr(invsyn(:,1),syn(:,1),'coeff')); % I corr inv

corr_syn_NCfr=max(xcorr(NCsyn(:,6),syn(:,6),'coeff')); % I correlation with raw amplitude regression
corr_syn_Cfr=max(xcorr(Csyn(:,6),syn(:,6),'coeff')); % I corr preconditioned regression
corr_syn_invfr=max(xcorr(invsyn(:,6),syn(:,6),'coeff')); % I corr inv


% used in frequency spectra plot 7
n_freq = 1/(2*ricker_sampling);
f_axis = (0:n_freq/(floor(nt/2)-1):n_freq);

%power spectra of the original synthetic, the noise spectrum, the
%conditioned ava spectrum, and the inverse spectrum, all with 6 columns of
%angle

%power_syn =   20*log10(bsxfun(@rdivide,abs(fft(syn_traces)),max(abs(fft(syn_traces)))));
%power_syn = power_syn(1:length(f_axis),:);
if coloured_noise==1
power_noise = 20*log10(bsxfun(@rdivide,abs(fft(syn_n_traces)),max(abs(fft(syn_n_traces)))));
power_noise = power_noise(1:length(f_axis),:);
end
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

%% results
% column 1 is NC, 2 is C, and 3 is inv
%after varying chosen variable, go to these matrices and copy into
%spreadsheet?  just save the param vectors and some key plots
%as i go, and reproduce the key results with more care using this script later
bpI = w_matrix(1:nt,1:nt)*I;
bpG = w_matrix(1:nt,1:nt)*G;
bpI = bpI*norm(I)/norm(bpI(subsample));
bpG = bpG*norm(G)/norm(bpG(subsample));

% res_corrI(zz,1)=corr_I_NC;
% res_corrI(zz,2)=corr_I_C;
% res_corrI(zz,3)=corr_I_inv;
% res_corrG(zz,1)=corr_G_NC;
% res_corrG(zz,2)=corr_G_C;
% res_corrG(zz,3)=corr_G_inv;
res_corrsynnr(zz,1)=corr_syn_NCnr;
res_corrsynnr(zz,2)=corr_syn_Cnr;
res_corrsynnr(zz,3)=corr_syn_invnr;
res_corrsynfr(zz,1)=corr_syn_NCfr;
res_corrsynfr(zz,2)=corr_syn_Cfr;
res_corrsynfr(zz,3)=corr_syn_invfr;


res_resi(zz,1)=norm(syn - (NCsyn*norm(syn)/norm(NCsyn)));
res_resi(zz,2)=norm(syn - (Csyn*norm(syn)/norm(Csyn)));
res_resi(zz,3)=norm(syn - (invsyn*norm(syn)/norm(invsyn)));

res_furthestIG(zz,1)=max(((scaledNCI - bpI).^2 + (scaledNCG - bpG).^2).^(0.5));
res_furthestIG(zz,2)=max(((scaledCI - bpI).^2 + (scaledCG - bpG).^2).^(0.5));
res_furthestIG(zz,3)=max(((scaledinvI - bpI).^2 + (scaledinvG - bpG).^2).^(0.5));
res_number_far(zz,1)=sum(((scaledNCI - bpI).^2 + (scaledNCG - bpG).^2).^(0.5)>0.1);
res_number_far(zz,2)=sum(((scaledCI - bpI).^2 + (scaledCG - bpG).^2).^(0.5)>0.1);
res_number_far(zz,3)=sum(((scaledinvI - bpI).^2 + (scaledinvG - bpG).^2).^(0.5)>0.1);

save('res.mat','res_corrsynnr','res_corrsynfr','res_resi','res_furthestIG','res_number_far');
end

if tvw==0
figure(1)
set(1,'Units','inches','Position',[0 0 23 16]);

elseif tvw==1
    figure(2)
    set(2,'Units','inches','Position',[0 0 23 16]);

end

axis_number_fontsize = 14;
axis_title_fontsize = 14;
title_fontsize = 16;
line_width = 2;
markersize=7;
subplot(2,2,1)
hold on
plot(res_corrsynnr(:,1),'^-.','Color',[1 0 0],'Markersize',markersize,'Linewidth',line_width)
plot(res_corrsynnr(:,2),'s-.','Color',[0 1 0],'Markersize',markersize,'Linewidth',line_width)
plot(res_corrsynnr(:,3),'o-.','Color',[0 0 1],'Markersize',markersize,'Linewidth',line_width)
plot(res_corrsynfr(:,1),'^-.','Color',[0.75 0 0],'Markersize',markersize,'Linewidth',line_width)
plot(res_corrsynfr(:,2),'s-.','Color',[0 0.75 0],'Markersize',markersize,'Linewidth',line_width)
plot(res_corrsynfr(:,3),'o-.','Color',[0 0 0.75],'Markersize',markersize,'Linewidth',line_width)
legend('NC near','C near','inv near','NC far','C far','inv far')
title('Correlation of real synthetic with estimated, nearest angle trace compared to furthest','FontSize',title_fontsize)
ylabel('Maximum correlation coefficient','FontSize',axis_title_fontsize)
xlabel('Test number','FontSize',axis_title_fontsize)

subplot(2,2,2)
hold on
plot(res_resi(:,1),'r^-.','Markersize',markersize,'Linewidth',line_width)
plot(res_resi(:,2),'gs-.','Markersize',markersize,'Linewidth',line_width)
plot(res_resi(:,3),'bo-.','Markersize',markersize,'Linewidth',line_width)
title('L_2 residual of synthetic compared to estimated','FontSize',title_fontsize)
xlabel('Test number','FontSize',axis_title_fontsize)
ylabel('Data-estimate L_2 residual','FontSize',axis_title_fontsize)
legend('NC','C','inv');

subplot(2,2,3)
hold on
plot(res_furthestIG(:,1),'r^-.','Markersize',markersize,'Linewidth',line_width)
plot(res_furthestIG(:,2),'gs-.','Markersize',markersize,'Linewidth',line_width)
plot(res_furthestIG(:,3),'bo-.','Markersize',markersize,'Linewidth',line_width)
title('Furthest estimated IG point from original location','FontSize',title_fontsize)
xlabel('Test number','FontSize',axis_title_fontsize)
ylabel('Furthest distance of estimated point compared to original','FontSize',axis_title_fontsize)
legend('NC','C','inv');

subplot(2,2,4)
hold on
plot(res_number_far(:,1),'r^-.','Markersize',markersize,'Linewidth',line_width)
plot(res_number_far(:,2),'gs-.','Markersize',markersize,'Linewidth',line_width)
plot(res_number_far(:,3),'bo-.','Markersize',markersize,'Linewidth',line_width)
title('Number of IG points at >0.001 euclidian distance from original','FontSize',title_fontsize)
xlabel('Test number','FontSize',axis_title_fontsize)
legend('NC','C','inv');
hold off
ylabel('Number','FontSize',axis_title_fontsize)
if tvw==0;
if save_on == 1
    set(gcf,'PaperPositionMode','auto')
    savefig('CNSN_stat_pw','eps','pdf','tiff','jpeg');
    saveas(figure(1),'CNSN_stat_pw.fig');
    hgsave('CNSN_stat','-v7.3');
    
end
else
    if save_on == 1
    set(gcf,'PaperPositionMode','auto')
    savefig('CNSN_tvw_pw','pdf','tiff');
    saveas(figure(2),'CNSN_tvw_pw.fig');
    hgsave('CNSN_tvw','-v7.3');
    
end
end

end

%% here is the scatter plots loop, choose which scatter to compare then make an if loop, scatter plots for one specific parameter test, 
% pick an interesting one, maybe beginning and end for a key test, compare and contrast
% if zz==
% if tvw==0
% figure(3)
% set(3,'Units','inches','Position',[0 0 23 16]);
% else
%     figure(4)
%     set(4,'Units','inches','Position',[0 0 23 16]);
%     
% end
% 
% bpI = w_matrix(1:nt,1:nt)*I;
% bpG = w_matrix(1:nt,1:nt)*G;
% bpI = bpI*norm(I)/norm(bpI(subsample));
% bpG = bpG*norm(G)/norm(bpG(subsample));
% 
% subplot(2,2,1)
% time = (1:1:nt);
% % time(IG==0) = 0;
% scatter(bpI(subsample),bpG(subsample),[],time(subsample),'filled','SizeData',50);
% axis([-0.3 0.3 -0.3 0.3]);
% colorbar
% t = colorbar('peer',gca);
% set(get(t,'ylabel'),'String', 'Sample Number','FontSize',axis_title_fontsize);
% title('True Intercept vs Gradient filtered to seismic bandwidth','FontSize',title_fontsize)
% xlabel('Intercept','FontSize',axis_title_fontsize)
% ylabel('Gradient','FontSize',axis_title_fontsize)
% set(gca,'FontSize',axis_number_fontsize)
% 
% subplot(2,2,2)
% % time = (1:1:nt);
% % time(NCIG==0) = 0;
% % scatter(NCI(subsample),NCG(subsample),[],time(subsample),'filled'); axis equal
% colours = ((scaledNCI - bpI).^2 + (scaledNCG - bpG).^2).^(0.5);
% scatter(scaledNCI(subsample),scaledNCG(subsample),[],colours(subsample),'filled','SizeData',50);
% axis([-0.3 0.3 -0.3 0.3]);
% caxis([0 0.25]);
% colorbar
% t = colorbar('peer',gca);
% set(get(t,'ylabel'),'String', 'Euclidean distance from true IG location','FontSize',axis_title_fontsize);
% title({'Intercept vs Gradient from raw amplitude regression';sprintf('Number of points >%.3f from true IG location = %d',0.1,sum(colours>0.1));sprintf('Furthest point from true IG location = %.3f',max(colours))},'FontSize',title_fontsize)
% xlabel('Intercept','FontSize',axis_title_fontsize)
% ylabel('Gradient','FontSize',axis_title_fontsize)
% set(gca,'FontSize',axis_number_fontsize)
% 
% subplot(2,2,3)
% % time = (1:1:nt);
% % time(CIG==0) = 0;
% % scatter(CI(subsample),CG(subsample),[],time(subsample),'filled'); axis equal
% colours = ((scaledCI - bpI).^2 + (scaledCG - bpG).^2).^(0.5);
% scatter(scaledCI(subsample),scaledCG(subsample),[],colours(subsample),'filled','SizeData',50);
% axis([-0.3 0.3 -0.3 0.3]);
% caxis([0 0.25]);
% colorbar
% t = colorbar('peer',gca);
% set(get(t,'ylabel'),'String', 'Euclidean distance from true IG location','FontSize',axis_title_fontsize);
% title({'Intercept vs Gradient from conditioned amplitude regression';sprintf('Number of points >%.3f from true IG location = %d',0.1,sum(colours>0.1));sprintf('Furthest point from true IG location = %.3f',max(colours))},'FontSize',title_fontsize)
% xlabel('Intercept','FontSize',axis_title_fontsize)
% ylabel('Gradient','FontSize',axis_title_fontsize)
% set(gca,'FontSize',axis_number_fontsize)
% 
% subplot(2,2,4)
% % time = (1:1:nt);
% % time(invIG==0) = 0;
% % scatter(invI(subsample),invG(subsample),[],time(subsample),'filled'); axis equal
% colours = ((scaledinvI - bpI).^2 + (scaledinvG - bpG).^2).^(0.5);
% scatter(scaledinvI(subsample),scaledinvG(subsample),[],colours(subsample),'filled','SizeData',50);
% axis([-0.3 0.3 -0.3 0.3]);
% caxis([0 0.25]);
% colorbar
% t = colorbar('peer',gca);
% set(get(t,'ylabel'),'String', 'Euclidean distance from true IG location','FontSize',axis_title_fontsize);
% title({'Intercept vs Gradient from inversion';sprintf('Number of points >%.3f from true IG location = %d',0.1,sum(colours>0.1));sprintf('Furthest point from true IG location = %.3f',max(colours))},'FontSize',title_fontsize)
% xlabel('Intercept','FontSize',axis_title_fontsize)
% ylabel('Gradient','FontSize',axis_title_fontsize)
% set(gca,'FontSize',axis_number_fontsize)
% if tvw==0
% if save_on == 1
%     set(gcf,'PaperPositionMode','auto')
%     savefig('crossplot_CNSN_stat','eps','pdf','tiff','jpeg');
%     hgsave('crossplot_CNSN_stat','-v7.3');
% end
% else
%     if save_on == 1
%     set(gcf,'PaperPositionMode','auto')
%     savefig('crossplot_CNSN_tvw','pdf','tiff','jpeg');
%     hgsave('crossplot_CNSN_tvw','-v7.3');
%     end
% end