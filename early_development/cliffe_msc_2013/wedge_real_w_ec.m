% This script is to simulate a wedge, so as to investigate the effects of
% offset dependent tuning
% this code makes a basic random synthetic, drops in a predetermined
% I and G at 2 specific locations, which vary on each loop, outputting a
% plot each time. only include randomisation steps in first loop, in all
% others, retain original parameters
close all
clear all

save_on = 1
if save_on ==1 % makes a new directory to save the pictures in
    cd('~/key plots');
end

%% make synthetic

nt = 1000;                                              % trace length in samples
sparsity = 0.7;                                         % proportion of non-zero amplitudes (0-1) 1 = fully sparse
max_amplitude = 0.15;
ricker_cf = [30 22 20 19 17 15];                        % central frequency of ricker wavelets to use to make synthetic traces (one per angle)
ricker_sampling = 0.004;                                % time sampling rate
background_dominace = 0.8;                              % proportion of the data that sits on the background IG trend (0-1) 1 = all background
chi = 21-(0:ricker_sampling:(nt-1)*ricker_sampling)';   % background trend angle (varies with time)
chinoise = 0;                                           % add noise to the chi angle (0-1) 0 = no noise, spread from the trends
IGnoise = 0;                                           % add noise to the IG data (0-1) 0 = no noise, general spread across cross plot
angles = [10 15 20 25 30 35];                           % angle of the traces
subsample = [1;zeros((0.004/ricker_sampling)-1,1)];
subsample = bsxfun(@times,ones(1,floor(nt/length(subsample))),subsample);
subsample = logical(subsample(:));
anomaly_multiplier = 1.5;                               % boost the off-trend amplitudes (>0)
tvw = 1; % switch for time varying wavelets (0 = off, 1 = on)

Q = 150;                                                % Q factor to be chosen                         
for i=1:length(angles)
tvw_Q(1,i) = exp(-(((pi()*ricker_cf(1,i)/cos(angles(1,i)/2*pi()*180)))/Q));     % the attenuation varies with frequency and incidence angle, with an assumption of t=1 and 0 degrees, t gets scaled larger for increasing offset              
end

tvw_att = tvw_Q;                                        % wavelet attenuation for time varying wavelets, = amplitude scaling for the final wavelet. need to assume a time/depth for the last trace to occur to do proper Q
IG_taper_length = 50;
coloured_noise = 1;    % switch for coloured noise (0 = off, 1 = on)
pre_white=1;
sn_ratio = [3,10,3];                                    % time varying signal to noise ratio (shallow,middle,deep)
cf_coloured_noise = [45 40 35 30 25 20];                % central frequency of angle-dependent coloured noise
for i=1:length(angles)
    tvw_att_n(1,i) = exp(-((pi()*cf_coloured_noise(1,i)/cos(angles(1,i)/2*pi()*180))/Q));      % Q for noise, as before. follows with the assumption that the noise is coherent, and attenuated. rather than random and uniform across angles
end

res_mov = 0;                                            % rnmo simulation not in this script

ref = synref(nt,sparsity,max_amplitude);                % make series of reflection coefficients, nt long
ref1=ref;
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

%code adjusts for shorter wedge matrix, so can put as many or as few
%intervals as you like
wedge(2,:)=[235,230,225,220,216,210,208,206,204,201];
   
  wedge(1,1:(length(wedge(2,:))))=200; 
  
for kk= 0:(length(wedge)-1)
    if kk>0;
        
       ref(wedge(2,kk),1)=ref1(wedge(2,1+kk),1);
    end
ref(wedge(1,1+kk),1)=0.25;
ref(wedge(2,1+kk),1)=0.25;
class(wedge(1,1+kk),1)=2;
class(wedge(2,1+kk),1)=2;

if kk==0;
R1=randn(nt,1);
R2=randn(nt,1);
R3=randn(nt,1);
R4=rand(2*nt,1);
else
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
        ar_tmp(:,ii)=w_tmp(1,1);   % this bit saves the original wavelet which is used to make the synthetic before time varying, so it can be compared to the wavelet estimated by the methods
        ar_tmp2(:,ii)=w_tmp(1,2);
        for j=1:size(ar_tmp{1,ii})
        real_w(j,ii)=ar_tmp{1,ii}(j,1);
        end
        for j=1:size(ar_tmp2{1,ii})
        real_w_att(j,ii)=ar_tmp2{1,ii}(j,1); %saves the attenuated original wavelet too
        end
        
        w_tmp = interp1([1;nt],[w_tmp{1}';w_tmp{2}'],(1:1:nt)');         % interpolates, between original wavelet and attenuated wavelet, fills huge matrix
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
    syn = w_matrix*[I;G];     %synthetic trace, wavelet*intercept and gradient

%Spectral Analysis of data without noise
for j=0:5            % these 2 loops take the synthetic trace from the 1 column the 6 angle gathers were in, and output into a 6 column matrix, so as to make the amplitude spectra
    for i=1:nt
syn_traces(i,j+1) = syn((i+j*nt),1);
    end
end

    syn_noise = noise_matrix*R4;      %noise matrix was previously wavelets of the central frequency given in the variables, now * random numbers, different ricker wavelets for angles, proportioned above, but same random numbers to proportion them
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

%plot(syn)
p_wedge{kk+1,1}=20*log10(bsxfun(@rdivide,abs(fft(syn_traces)),max(abs(fft(syn_traces))))); %power spectra of each angle gather, will be looped over

%% solve non-conditioned AVA

NCava = [ones(size(angles')) sin(angles'*pi/180).*sin(angles'*pi/180)]\syn';

NCI = NCava(1,:)';
NCG = NCava(2,:)';
NCIG = NCI.*NCG;
NCIG(NCIG<0) = 0;

%% solve conditioned AVA
for ii=1:length(angles)

w_est(:,ii)=full(w_matrix(((nt/2)-30+nt*(ii-1)):((nt/2)+30+(nt*(ii-1))),nt/2));%feeding the average wavelet into the matching filter production, position in w_matrix irrelevent when tvw=0, but picks middle interpolated wavelet (.: average?) when tvw=1

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

Cava = [ones(size(angles')) sin(angles'*pi/180).*sin(angles'*pi/180)]\Csyn'; % inversion

CI = Cava(1,:)';
CG = Cava(2,:)';
CIG = CI.*CG;
CIG(CIG<0) = 0;

%% Inverted AVA
wsmooth = 0;                                    % variable, controls smoothness, but reduces how good the fit is, never do this on real data prior to checking the unsmoothed data
n_wavelets = 3;

smooth = spdiags([-wsmooth*ones(2*nt,1) 2*wsmooth*ones(2*nt,1) -wsmooth*ones(2*nt,1)],[-1 0 1],2*nt,2*nt);
    IGmatrix = [w_matrix; smooth];


if pre_white==1;
pw = spdiags(ones(nt,2),[0,nt],nt,2*nt);
tmp = sind(angles).*sind(angles);
tmp = repmat(tmp,nt,1);
tmp = tmp(:);
pw = repmat(pw,length(angles),1);
pw(:,1001:end) = bsxfun(@times,pw(:,1001:end),(tmp*0.5));
pw = [pw;zeros(2*nt,2*nt)];
    
    IGmatrix = IGmatrix+pw;
else
end
data = [syn(:);zeros(2*nt,1)];
model = [NCI;NCG];

invava = lsqr(IGmatrix,data,1e-8,2000,[],[],model);  %AVA inversion, uses unconditioned AVA inversion as starting model, and synthetic as input data

invI = invava(1:nt); 
invG = invava(nt+1:end);

invIG = invI.*invG;
invIG(invIG<0) = 0;
%% syn plots

axis_number_fontsize = 14;
axis_title_fontsize = 14;
title_fontsize = 16;
line_width = 1;
% 
% figure(1)
% set(1,'Units','inches','Position',[0 0 28 10]);
% 
% 
% subplot(length(wedge),1,1+kk)
% plot(syn,'LineWidth',line_width)
% xlabel('Sample number','FontSize',axis_title_fontsize)
% ylabel('Amplitude','FontSize',axis_title_fontsize)
% set(gca,'FontSize',axis_number_fontsize)
% axis([180 260 -0.4 0.4])
% 
% if save_on == 1
%     set(gcf,'PaperPositionMode','auto');
%     savefig('wedge_synthetic','eps','pdf','tiff','jpeg');
%     hgsave('wedge_synthetic','-v7.3');
% end


%% residual plots

scaledNCI = NCI*norm(I)/norm(NCI(subsample));
scaledNCG = NCG*norm(G)/norm(NCG(subsample));
scaledCI = CI*norm(I)/norm(CI(subsample));
scaledCG = CG*norm(G)/norm(CG(subsample));
scaledinvI = invI*norm(I)/norm(invI(subsample));
scaledinvG = invG*norm(G)/norm(invG(subsample));

NCsyn=reshape(w_matrix(1:nt*length(angles),:)*[NCI;NCG],nt,[]);
Csyn= reshape(w_matrix(1:nt*length(angles),:)*[CI;CG],nt,[]);
invsyn= reshape(w_matrix(1:nt*length(angles),:)*[invI;invG],nt,[]);

Vplot(kk+1)=norm(syn-(NCsyn*norm(syn)/norm(NCsyn)));
Wplot(kk+1)=norm(syn-(Csyn*norm(syn)/norm(Csyn)));
Xplot(kk+1)=norm(syn-(invsyn*norm(syn)/norm(invsyn)));

%% figure 3 & 4

Iwedge(kk+1,1)=I(wedge(2,kk+1),1);
NCIwedge(kk+1,1)=NCI(wedge(2,kk+1),1);
CIwedge(kk+1,1)=CI(wedge(2,kk+1),1);
invIwedge(kk+1,1)=invI(wedge(2,kk+1),1);

Gwedge(kk+1,1)=G(wedge(2,kk+1),1);
NCGwedge(kk+1,1)=NCG(wedge(2,kk+1),1);
CGwedge(kk+1,1)=CG(wedge(2,kk+1),1);
invGwedge(kk+1,1)=invG(wedge(2,kk+1),1);

end

axis_number_fontsize = 14;
axis_title_fontsize = 14;
title_fontsize = 16;
line_width = 2;
markersize=7;
figure(2)
    hold
plot((wedge(2,:)-wedge(1,:)),Vplot,'r^-.','Linewidth',line_width,'Markersize',markersize)
title('Plot of resiudals with decreasing bed thickness','FontSize',title_fontsize)
xlabel('Wedge thickness measured in traces','FontSize',axis_title_fontsize)
ylabel('residual of IG estimate','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)
plot((wedge(2,:)-wedge(1,:)),Wplot,'gs-.','Linewidth',line_width,'Markersize',markersize)
plot((wedge(2,:)-wedge(1,:)),Xplot,'bo-.','Linewidth',line_width,'Markersize',markersize)
LEGEND('non-conditioned estimate','conditioned data estimate','DIGI')

if save_on == 1
    set(gcf,'PaperPositionMode','auto')
    savefig('wedge_rw_residuals','pdf','tiff');
    hgsave('wedge_rw residuals','-v7.3');
    saveas(figure(2),'ang dep phase tvw.fig');

end


figure(3)
hold
title('Plot of I against increasing wedge thickness')
plot((wedge(2,:)-wedge(1,:)),Iwedge,'k','Linewidth',line_width,'Markersize',markersize)
plot((wedge(2,:)-wedge(1,:)),NCIwedge,'r^-.','Linewidth',line_width,'Markersize',markersize)
plot((wedge(2,:)-wedge(1,:)),CIwedge,'gs-.','Linewidth',line_width,'Markersize',markersize)
plot((wedge(2,:)-wedge(1,:)),invIwedge,'bo-.','Linewidth',line_width,'Markersize',markersize)
xlabel('Wedge thickness measured in traces','FontSize',axis_title_fontsize)
ylabel('Value of I','FontSize',axis_title_fontsize)
LEGEND('real I data','non-conditioned estimate','conditioned data estimate','DIGI')
    if save_on == 1
    set(gcf,'PaperPositionMode','auto')
    savefig('wedge_rw I','pdf','tiff');
    saveas(figure(3),'wedge_rw I.fig');
    hgsave('wedge_rw I','-v7.3');
    end
figure(4)
hold
title('Plot of G against increasing wedge thickness')
plot((wedge(2,:)-wedge(1,:)),Gwedge,'k','Linewidth',line_width,'Markersize',markersize)
plot((wedge(2,:)-wedge(1,:)),NCGwedge,'r^-.','Linewidth',line_width,'Markersize',markersize)
plot((wedge(2,:)-wedge(1,:)),CGwedge,'gs-.','Linewidth',line_width,'Markersize',markersize)
plot((wedge(2,:)-wedge(1,:)),invGwedge,'bo-.','Linewidth',line_width,'Markersize',markersize)
LEGEND('real G data','non-conditioned estimate','conditioned data estimate','DIGI')
xlabel('Wedge thickness measured in traces','FontSize',axis_title_fontsize)
ylabel('Value of G','FontSize',axis_title_fontsize)
    if save_on == 1
    set(gcf,'PaperPositionMode','auto')
    savefig('wedge_rw G','pdf','tiff');
    saveas(figure(4),'wedge_rw G.fig');
    hgsave('wedge_rw G','-v7.3');
    end