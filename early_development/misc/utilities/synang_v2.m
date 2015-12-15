close all
clear all

save_on = 0; % 0 saves no pictures, 1 saves all pictures

if save_on ==1 % makes a new directory to save the pictures in
    test_id = round(10000*rand);
    system(sprintf('mkdir digi_test_%d',test_id));
    cd(sprintf('digi_test_%d',test_id));
end
%% make synthetic

nt = 1000;                                              % trace length in samples
sparsity = 0.7;                                         % proportion of non-zero amplitudes (0-1) 1 = fully sparse
max_amplitude = 0.15;
ricker_cf = [27 22 20 19 17 15];                        % central frequency of ricker wavelets to use to make synthetic traces (one per angle)
ricker_sampling = 0.004;                                % time sampling rate
background_dominace = 0.8;                              % proportion of the data that sits on the background IG trend (0-1) 1 = all background
chi = 21-(0:ricker_sampling:(nt-1)*ricker_sampling)';   % background trend angle (varies with time)
chinoise = 0;                                           % add noise to the chi angle (0-1) 0 = no noise
IGnoise = 0.8;                                          % add noise to the IG data (0-1) 0 = no noise
angles = [10 15 20 25 30 35];                           % angle of the traces
subsample = [1;zeros((0.004/ricker_sampling)-1,1)];
subsample = bsxfun(@times,ones(1,floor(nt/length(subsample))),subsample);
subsample = logical(subsample(:));
anomaly_multiplier = 1.5;                               % boost the off-trend amplitudes (>0)
tvw = 0;                                                % switch for time varying wavelets (0 = off, 1 = on)
tvw_att = 0.4;                                          % wavelet attenuation for time varying wavelets
IG_taper_length = 50;
coloured_noise = 0;                                     % switch for coloured noise (0 = off, 1 = on)
sn_ratio = [3,10,3];                                    % time varying signal to noise ratio (shallow,middle,deep)
cf_coloured_noise = [45 40 35 30 25 20];                % central frequency of angle-dependent coloured noise

ref = synref(nt,sparsity,max_amplitude);                % make reflection coefficients
class = nan(nt,1);
random = rand(nt,1);
class(ref~=0) = random(ref~=0);
class(class>=(1-background_dominace)) = 1;
class(class<(1-background_dominace)) = 2;
class_2_flag = 0;
for ii = 1:nt
    if class(ii) == 2
        class_2_flag = 1;
        ref(ii) = -anomaly_multiplier*abs(ref(ii));
    end
    if class(ii) == 1
        if class_2_flag == 1
            class(ii) = 3;
            ref(ii) = anomaly_multiplier*abs(ref(ii));
            class_2_flag = 0;
        end
    end
end

I = ref;
G = zeros(nt,1);
chi_array = chi+(chinoise-2*chinoise*rand(nt,1));
G(class==1) = -I(class==1)./tan(chi_array(class==1)*pi/180);
G(class==2) = -I(class==2)./tan(90+chi_array(class==2)*pi/180);
G(class==3) = -I(class==3)./tan(90+chi_array(class==3)*pi/180);

Gnoise = IGnoise*std(G)*randn(nt,1);
G(~isnan(class)) = G(~isnan(class))+Gnoise(~isnan(class));

Inoise = IGnoise*std(I)*randn(nt,1);
I(~isnan(class)) = I(~isnan(class))+Inoise(~isnan(class));

I = bsxfun(@times,[fliplr(0.5+0.5*cos(0:pi/IG_taper_length:pi)),ones(1,nt-2*IG_taper_length-2),0.5+0.5*cos(0:pi/IG_taper_length:pi)]',I);
G = bsxfun(@times,[fliplr(0.5+0.5*cos(0:pi/IG_taper_length:pi)),ones(1,nt-2*IG_taper_length-2),0.5+0.5*cos(0:pi/IG_taper_length:pi)]',G);

eref = bsxfun(@plus,I,bsxfun(@times,G,sin(angles*pi/180).*sin(angles*pi/180)));

if tvw == 1
    for ii = 1:length(angles)
        w_tmp = {ricker(ricker_cf(ii),ricker_sampling),tvw_att*ricker(tvw_att*ricker_cf(ii),ricker_sampling)};
        w_tmp{1} = circshift([w_tmp{1};zeros(length(w_tmp{2})-length(w_tmp{1}),1)],(length(w_tmp{2})-length(w_tmp{1}))/2);
        w_tmp = interp1([1;nt],[w_tmp{1}';w_tmp{2}'],(1:1:nt)');
        if ii == 1
            w_matrix = [spdiags(w_tmp,-floor(size(w_tmp,2)/2):floor(size(w_tmp,2)/2),nt,nt), spdiags(bsxfun(@times,sind(angles(ii))*sind(angles(ii)),w_tmp),-floor(size(w_tmp,2)/2):floor(size(w_tmp,2)/2),nt,nt)];
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
            w_tmp = {ricker(cf_coloured_noise(ii),ricker_sampling),tvw_att*ricker(tvw_att*cf_coloured_noise(ii),ricker_sampling)};
            w_tmp{1} = circshift([w_tmp{1};zeros(length(w_tmp{2})-length(w_tmp{1}),1)],(length(w_tmp{2})-length(w_tmp{1}))/2);
            w_tmp = interp1([1;nt],[w_tmp{1}';w_tmp{2}'],(1:1:nt)');
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
    syn = w_matrix*[I;G];
    syn_noise = noise_matrix*rand(2*nt,1);
    syn_noise = (syn_noise./repmat(sqrt(interp1([1 round(nt/2) nt],sn_ratio,(1:1:nt)','linear')),length(angles),1))*(std(syn)/std(syn_noise));
    syn = reshape(syn+syn_noise,nt,[]);
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

ft_taper_length = 10; % assymetric (tapers only the high frequency energy to zero at Nyquist frequency)
wavelet_taper_length = 30; % symmetric about zero (outside +/- wavelet taper length about zero, amplitude is set to zero, while inside there is a cosine taper)
wavelet_length = 151;
half_wavelet_length = floor(wavelet_length/2);
hnt = floor(nt/2);
ft_taper = [ones(hnt-ft_taper_length,1); ((1+cos((0:1/(ft_taper_length-1):1)*pi)')/2)];
ft_taper = [ft_taper; zeros(nt-length(ft_taper),1)];

ft_traces_tmp = abs(fft(syn));
ft_traces_tmp = bsxfun(@times,ft_traces_tmp,ft_taper); % ensures zero energy at Nyquist frequency
w_est = circshift(ifft(ft_traces_tmp,'symmetric'),hnt);
[~, peak_index] = max(w_est);
peak_index = max(peak_index);
w_est = w_est(peak_index-half_wavelet_length:peak_index+half_wavelet_length,:);
wavelet_taper = [((1+cos((-1:1/(wavelet_taper_length-1):0)*pi)')/2); ones(13,1); ((1+cos((0:1/(wavelet_taper_length-1):1)*pi)')/2)];
wavelet_taper = [zeros((size(w_est,1)-length(wavelet_taper))/2,1); wavelet_taper; zeros((size(w_est,1)-length(wavelet_taper))/2,1)];
w_est = bsxfun(@times,w_est,wavelet_taper); % ensures zero amplitude at terminations
w_est(isnan(w_est)) = 0;

for ii = 1:length(angles)
    filter(:,ii) = match_call(w_est(:,round((length(angles))/2)),w_est(:,ii),-6);
end

for ii = 1:length(angles)
    Csyn(:,ii) = conv(syn(:,ii),filter(:,ii),'same');
    fscalar = norm(syn(:,ii))/norm(Csyn(:,ii));
    filter(:,ii) = filter(:,ii)*fscalar;
    Csyn(:,ii) = conv(syn(:,ii),filter(:,ii),'same');
end

Cava = [ones(size(angles')) sin(angles'*pi/180).*sin(angles'*pi/180)]\Csyn';

CI = Cava(1,:)';
CG = Cava(2,:)';
CIG = CI.*CG;
CIG(CIG<0) = 0;

%% Inverted AVA
wsmooth = 0;
n_wavelets = 3;

smooth = spdiags([-wsmooth*ones(2*nt,1) 2*wsmooth*ones(2*nt,1) -wsmooth*ones(2*nt,1)],[-1 0 1],2*nt,2*nt);

if tvw == 1
    wavelet_step = floor(nt/(n_wavelets+1));
    wavelet_grid = (wavelet_step:wavelet_step:n_wavelets*wavelet_step);
    ft_taper_length = 10; % assymetric (tapers only the high frequency energy to zero at Nyquist frequency)
    wavelet_taper_length = 30; % symmetric about zero (outside +/- wavelet taper length about zero, amplitude is set to zero, while inside there is a cosine taper)
    wavelet_length = 301;
    half_wavelet_length = floor(wavelet_length/2);
    hnt = wavelet_step;
    ft_taper = [ones(hnt-ft_taper_length,1); ((1+cos((0:1/(ft_taper_length-1):1)*pi)')/2)];
    ft_taper = [ft_taper; zeros(2*wavelet_step-length(ft_taper),1)];
    
    for jj = 1:n_wavelets
        ft_traces_tmp = abs(fft(syn(1+(jj-1)*wavelet_step:(jj+1)*wavelet_step,:)));
        ft_traces_tmp = bsxfun(@times,ft_traces_tmp,ft_taper); % ensures zero energy at Nyquist frequency
        w_est = circshift(ifft(ft_traces_tmp,'symmetric'),hnt);
        [~, peak_index] = max(w_est);
        peak_index = max(peak_index);
        w_est = w_est(peak_index-half_wavelet_length:peak_index+half_wavelet_length,:);
        wavelet_taper = [((1+cos((-1:1/(wavelet_taper_length-1):0)*pi)')/2); ones(13,1); ((1+cos((0:1/(wavelet_taper_length-1):1)*pi)')/2)];
        wavelet_taper = [zeros((size(w_est,1)-length(wavelet_taper))/2,1); wavelet_taper; zeros((size(w_est,1)-length(wavelet_taper))/2,1)];
        w_est = bsxfun(@times,w_est,wavelet_taper); % ensures zero amplitude at terminations
        w_est(isnan(w_est)) = 0;
        normalise = sum(sqrt(w_est.^2));
        w_est = normalise(1)*bsxfun(@rdivide,w_est,normalise);
        w_est = w_est(:);
        w_est_tmp(:,jj) = w_est;
    end
    
    w_est = w_est_tmp;
    
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

invava = lsqr(IGmatrix,data,1e-3,500,[],[],model);

invI = invava(1:nt);
invG = invava(nt+1:end);

invIG = invI.*invG;
invIG(invIG<0) = 0;

%% plots

axis_number_fontsize = 16;
axis_title_fontsize = 16;
title_fontsize = 16;
line_width = 2;


%% syn plots

scaledNCI = NCI*norm(I)/norm(NCI(subsample));
scaledNCG = NCG*norm(G)/norm(NCG(subsample));
scaledCI = CI*norm(I)/norm(CI(subsample));
scaledCG = CG*norm(G)/norm(CG(subsample));
scaledinvI = invI*norm(I)/norm(invI(subsample));
scaledinvG = invG*norm(G)/norm(invG(subsample));

figure(1)

set(1,'Units','inches','Position',[0 0 28 10]);

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

%% results plot

figure(2)

set(2,'Units','inches','Position',[0 0 28 10]);

subplot(3,1,1)
stem(I,'marker','none','color',[150/255 150/255 150/255],'LineWidth',line_width)
hold all
plot(scaledNCI,'LineWidth',line_width)
hold off
axis tight
title('Intercept from raw amplitude regression (truth in grey)','FontSize',title_fontsize)
xlabel('Sample number','FontSize',axis_title_fontsize)
ylabel('Amplitude','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

subplot(3,1,2)
stem(I,'marker','none','color',[150/255 150/255 150/255],'LineWidth',line_width)
hold all
plot(scaledCI,'LineWidth',line_width)
hold off
axis tight
title('Intercept from conditioned amplitude regression (truth in grey)','FontSize',title_fontsize)
xlabel('Sample number','FontSize',axis_title_fontsize)
ylabel('Amplitude','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

subplot(3,1,3)
stem(I,'marker','none','color',[150/255 150/255 150/255],'LineWidth',line_width)
hold all
plot(scaledinvI,'LineWidth',line_width)
hold off
axis tight
title('Intercept from inversion (truth in grey)','FontSize',title_fontsize)
xlabel('Sample number','FontSize',axis_title_fontsize)
ylabel('Amplitude','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

if save_on == 1
    set(gcf,'PaperPositionMode','auto')
    savefig('intercept','eps','pdf','tiff','jpeg');
    hgsave('intercept','-v7.3');
end

figure(3)

set(3,'Units','inches','Position',[0 0 28 10]);

subplot(3,1,1)
stem(G,'marker','none','color',[150/255 150/255 150/255],'LineWidth',line_width)
hold all
plot(scaledNCG,'LineWidth',line_width)
hold off
axis tight
title('Gradient from raw amplitude regression (truth in grey)','FontSize',title_fontsize)
xlabel('Sample number','FontSize',axis_title_fontsize)
ylabel('Amplitude','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

subplot(3,1,2)
stem(G,'marker','none','color',[150/255 150/255 150/255],'LineWidth',line_width)
hold all
plot(scaledCG,'LineWidth',line_width)
hold off
axis tight
title('Gradient from conditioned amplitude regression (truth in grey)','FontSize',title_fontsize)
xlabel('Sample number','FontSize',axis_title_fontsize)
ylabel('Amplitude','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

subplot(3,1,3)
stem(G,'marker','none','color',[150/255 150/255 150/255],'LineWidth',line_width)
hold all
plot(scaledinvG,'LineWidth',line_width)
hold off
axis tight
title('Gradient from inversion (truth in grey)','FontSize',title_fontsize)
xlabel('Sample number','FontSize',axis_title_fontsize)
ylabel('Amplitude','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

if save_on == 1
    set(gcf,'PaperPositionMode','auto')
    savefig('gradient','eps','pdf','tiff','jpeg');
    hgsave('gradient','-v7.3');
end


%% scatter plots

figure(4)

set(4,'Units','inches','Position',[0 0 23 16]);

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

%% eer plots

figure(5)

set(5,'Units','inches','Position',[0 0 28 10]);

subplot(3,1,1)
stem(I.*cosd(chi) + G.*sind(chi),'marker','none','color',[150/255 150/255 150/255],'LineWidth',line_width)
hold all
plot(scaledNCI.*cosd(chi) + scaledNCG.*sind(chi),'LineWidth',line_width)
hold off
axis tight
title('Minimum energy EER projection using I and G from raw amplitude regression (truth in grey)','FontSize',title_fontsize)
xlabel('Sample Number','FontSize',axis_title_fontsize)
ylabel('Amplitude','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

subplot(3,1,2)
stem(I.*cosd(chi) + G.*sind(chi),'marker','none','color',[150/255 150/255 150/255],'LineWidth',line_width)
hold all
plot(scaledCI.*cosd(chi) + scaledCG.*sind(chi),'LineWidth',line_width)
hold off
axis tight
title('Minimum energy EER projection using I and G from conditioned amplitude regression (truth in grey)','FontSize',title_fontsize)
xlabel('Sample Number','FontSize',axis_title_fontsize)
ylabel('Amplitude','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

subplot(3,1,3)
stem(I.*cosd(chi) + G.*sind(chi),'marker','none','color',[150/255 150/255 150/255],'LineWidth',line_width)
hold all
plot(scaledinvI.*cosd(chi) + scaledinvG.*sind(chi),'LineWidth',line_width)
hold off
axis tight
title('Minimum energy EER projection using I and G from inversion (truth in grey)','FontSize',title_fontsize)
xlabel('Sample Number','FontSize',axis_title_fontsize)
ylabel('Amplitude','FontSize',axis_title_fontsize)
set(gca,'FontSize',axis_number_fontsize)

if save_on == 1
    set(gcf,'PaperPositionMode','auto')
    savefig('eerme','eps','pdf','tiff','jpeg');
    hgsave('eerme','-v7.3');
end


%% residual plots

figure(6)

set(6,'Units','inches','Position',[0 0 23 16]);

ncolours = 64;
colours = [[ones(ncolours/2,1);(1:-1/((ncolours/2)-1):0)'] [(0:1/((ncolours/2)-2):1)';1;1;(1:-1/((ncolours/2)-2):0)'] [(0:1/((ncolours/2)-1):1)';ones(ncolours/2,1)]];

NCsyn = reshape(w_matrix(1:nt*length(angles),:)*[NCI;NCG],nt,[]);
Csyn = reshape(w_matrix(1:nt*length(angles),:)*[CI;CG],nt,[]);
invsyn = reshape(w_matrix(1:nt*length(angles),:)*[invI;invG],nt,[]);

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