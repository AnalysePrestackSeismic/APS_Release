 run synang_param_real_w_ec

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
% Icw=conv(I,w_est_ava(:,1),'same');
% Gcw=conv(G,w_est_ava(:,1),'same');
% Iinvw=conv(I,w_est_tr(:,1),'same');
% Ginvw=conv(G,w_est_tr(:,1),'same');
% 
% corr_I_NC=max(xcorr(Incw,scaledNCI,'coeff')); % I correlation with raw amplitude regression
% corr_G_NC=max(xcorr(Gncw,scaledNCG,'coeff')); % G correlation with raw regression
% corr_I_C=max(xcorr(Icw,scaledCI,'coeff')); % I corr preconditioned regression
% corr_G_C=max(xcorr(Gcw,scaledCG,'coeff')); % G corr preconditioned regression
% corr_I_inv=max(xcorr(Iinvw,scaledinvI,'coeff')); % I corr inv
% corr_G_inv=max(xcorr(Ginvw,scaledinvG,'coeff')); % G corr inv

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
% power_w_est_C = 20*log10(bsxfun(@rdivide,abs(fft(w_est_ava)),max(abs(fft(w_est_ava)))));
% f_c= (0:n_freq/((wavelet_length_C/2)-1):n_freq);
% power_w_est_inv = 20*log10(bsxfun(@rdivide,abs(fft(w_est_tr)),max(abs(fft(w_est_tr)))));
% f_inv= (0:n_freq/((wavelet_length_inv/2)-1):n_freq);



%% plots

axis_number_fontsize = 14;
axis_title_fontsize = 14;
title_fontsize = 16;
line_width = 2;


%% syn plots

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
% 
% if save_on == 1
%     set(gcf,'PaperPositionMode','auto');
%     savefig('synthetic','eps','pdf','tiff','jpeg');
%     hgsave('synthetic','-v7.3');
% end
% 
% % results plot
% 
% figure(2)
% 
% set(2,'Units','inches','Position',[0 0 28 10]);
% 
% subplot(3,1,1)
% stem(I,'marker','none','color',[150/255 150/255 150/255],'LineWidth',line_width)
% hold all
% plot(scaledNCI,'LineWidth',line_width)
% hold off
% axis tight
% title(sprintf('Intercept from raw amplitude regression (truth in grey) %1.3g',corr_I_NC), 'FontSize',title_fontsize)
% xlabel('Sample number','FontSize',axis_title_fontsize)
% ylabel('Amplitude','FontSize',axis_title_fontsize)
% set(gca,'FontSize',axis_number_fontsize)
% 
% subplot(3,1,2)
% stem(I,'marker','none','color',[150/255 150/255 150/255],'LineWidth',line_width)
% hold all
% plot(scaledCI,'LineWidth',line_width)
% hold off
% axis tight
% title(sprintf('Intercept from conditioned amplitude regression (truth in grey) %1.3g',corr_I_C),'FontSize',title_fontsize)
% xlabel('Sample number','FontSize',axis_title_fontsize)
% ylabel('Amplitude','FontSize',axis_title_fontsize)
% set(gca,'FontSize',axis_number_fontsize)
% 
% subplot(3,1,3)
% stem(I,'marker','none','color',[150/255 150/255 150/255],'LineWidth',line_width)
% hold all
% plot(scaledinvI,'LineWidth',line_width)
% hold off
% axis tight
% title(sprintf('Intercept from inversion (truth in grey) =%1.3g',corr_I_inv),'FontSize',title_fontsize)
% xlabel('Sample number','FontSize',axis_title_fontsize)
% ylabel('Amplitude','FontSize',axis_title_fontsize)
% set(gca,'FontSize',axis_number_fontsize)
% 
% if save_on == 1
%     set(gcf,'PaperPositionMode','auto')
%     savefig('intercept','eps','pdf','tiff','jpeg');
%     hgsave('intercept','-v7.3');
% end
% 
% figure(3)
% 
% set(3,'Units','inches','Position',[0 0 28 10]);
% 
% subplot(3,1,1)
% stem(G,'marker','none','color',[150/255 150/255 150/255],'LineWidth',line_width)
% hold all
% plot(scaledNCG,'LineWidth',line_width)
% hold off
% axis tight
% title(sprintf('Gradient from raw amplitude regression (truth in grey) %1.3g', corr_G_NC),'FontSize',title_fontsize)
% xlabel('Sample number','FontSize',axis_title_fontsize)
% ylabel('Amplitude','FontSize',axis_title_fontsize)
% set(gca,'FontSize',axis_number_fontsize)
% 
% subplot(3,1,2)
% stem(G,'marker','none','color',[150/255 150/255 150/255],'LineWidth',line_width)
% hold all
% plot(scaledCG,'LineWidth',line_width)
% hold off
% axis tight
% title(sprintf('Gradient from conditioned amplitude regression (truth in grey) %1.3g',corr_G_C),'FontSize',title_fontsize)
% xlabel('Sample number','FontSize',axis_title_fontsize)
% ylabel('Amplitude','FontSize',axis_title_fontsize)
% set(gca,'FontSize',axis_number_fontsize)
% 
% subplot(3,1,3)
% stem(G,'marker','none','color',[150/255 150/255 150/255],'LineWidth',line_width)
% hold all
% plot(scaledinvG,'LineWidth',line_width)
% hold off
% axis tight
% title(sprintf('Gradient from inversion (truth in grey) %1.3g',corr_G_inv),'FontSize',title_fontsize)
% xlabel('Sample number','FontSize',axis_title_fontsize)
% ylabel('Amplitude','FontSize',axis_title_fontsize)
% set(gca,'FontSize',axis_number_fontsize)
% 
% if save_on == 1
%     set(gcf,'PaperPositionMode','auto')
%     savefig('gradient','eps','pdf','tiff','jpeg');
%     hgsave('gradient','-v7.3');
% end
% 
% 
% scatter plots
% 
% figure(4)
% 
% set(4,'Units','inches','Position',[0 0 23 16]);
% 
% bpI = w_matrix(1:nt,1:nt)*I;
% bpG = w_matrix(1:nt,1:nt)*G;
% bpI = bpI*norm(I)/norm(bpI(subsample));
% bpG = bpG*norm(G)/norm(bpG(subsample));
% 
% subplot(2,2,1)
% time = (1:1:nt);
% time(IG==0) = 0;
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
% time = (1:1:nt);
% time(NCIG==0) = 0;
% scatter(NCI(subsample),NCG(subsample),[],time(subsample),'filled'); axis equal
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
% time = (1:1:nt);
% time(CIG==0) = 0;
% scatter(CI(subsample),CG(subsample),[],time(subsample),'filled'); axis equal
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
% time = (1:1:nt);
% time(invIG==0) = 0;
% scatter(invI(subsample),invG(subsample),[],time(subsample),'filled'); axis equal
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
% 
% if save_on == 1
%     set(gcf,'PaperPositionMode','auto')
%     savefig('crossplot','eps','pdf','tiff','jpeg');
%     hgsave('crossplot','-v7.3');
% end
% 
% % eer plots
% 
% figure(5)
% 
% set(5,'Units','inches','Position',[0 0 28 10]);
% 
% subplot(3,1,1)
% stem(I.*cosd(chi) + G.*sind(chi),'marker','none','color',[150/255 150/255 150/255],'LineWidth',line_width)
% hold all
% plot(scaledNCI.*cosd(chi) + scaledNCG.*sind(chi),'LineWidth',line_width)
% hold off
% axis tight
% title('Minimum energy EER projection using I and G from raw amplitude regression (truth in grey)','FontSize',title_fontsize)
% xlabel('Sample Number','FontSize',axis_title_fontsize)
% ylabel('Amplitude','FontSize',axis_title_fontsize)
% set(gca,'FontSize',axis_number_fontsize)
% 
% subplot(3,1,2)
% stem(I.*cosd(chi) + G.*sind(chi),'marker','none','color',[150/255 150/255 150/255],'LineWidth',line_width)
% hold all
% plot(scaledCI.*cosd(chi) + scaledCG.*sind(chi),'LineWidth',line_width)
% hold off
% axis tight
% title('Minimum energy EER projection using I and G from conditioned amplitude regression (truth in grey)','FontSize',title_fontsize)
% xlabel('Sample Number','FontSize',axis_title_fontsize)
% ylabel('Amplitude','FontSize',axis_title_fontsize)
% set(gca,'FontSize',axis_number_fontsize)
% 
% subplot(3,1,3)
% stem(I.*cosd(chi) + G.*sind(chi),'marker','none','color',[150/255 150/255 150/255],'LineWidth',line_width)
% hold all
% plot(scaledinvI.*cosd(chi) + scaledinvG.*sind(chi),'LineWidth',line_width)
% hold off
% axis tight
% title('Minimum energy EER projection using I and G from inversion (truth in grey)','FontSize',title_fontsize)
% xlabel('Sample Number','FontSize',axis_title_fontsize)
% ylabel('Amplitude','FontSize',axis_title_fontsize)
% set(gca,'FontSize',axis_number_fontsize)
% 
% if save_on == 1
%     set(gcf,'PaperPositionMode','auto')
%     savefig('eerme','pdf','tiff');
%     saveas(figure(5),'eerme.fig');
%     hgsave('eerme','-v7.3');
% end
% 
% 
% residual plots

figure(6)

set(6,'Units','inches','Position',[0 0 23 16]);

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
    savefig('gathers_resmovmax','eps','pdf','tiff','jpeg');
    hgsave('gathers_resmovmax','-v7.3');
saveas(figure(6),'gathers_resmovmax.fig');
    save('workspace.mat','-v7.3');
end

%% Frequency Spectra plots
%these plots compare the angle stacks of the original data amplitude
%spectra with noise, and the noisy conditioned data amplitude spectra
% 
% axis_number_fontsize = 14;
% axis_title_fontsize = 16;
% title_fontsize = 17;
% line_width = 0.5;
% line_width2=2;
% marker_size=8;
% % 
% figure (7)
% set(7,'Units','inches','Position',[0 0 23 16]);
% for j=1:6
%      corr_noise(j,1)=max(xcorr(power_syn(:,j),power_noise(:,j),'coeff'));  % correlations for each angle gather, noise and original
%      corr_Csyn(j,1)=max(xcorr(power_syn(:,j),power_Csyn(:,j),'coeff'));    % conditioned synthetic and original
%      corr_inv(j,1)=max(xcorr(power_syn(:,j),power_invsyn(:,j),'coeff'));   % derived synthetic using estimated I and G and original correct wavelets, compared with original
%      corr_inv_n(j,1)=max(xcorr(power_invsyn(:,j),power_noise(:,j),'coeff')); %correlation of inversion synthetic with noise
%      corr_Csyn_n(j,1)=max(xcorr(power_Csyn(:,j),power_noise(:,j),'coeff')); % correlation of conditioned synthetic with noise
% end
%  subplot('position',[0.1 0.750 0.8 0.12])
%  hold on
%  plot(f_axis,power_noise(:,1),'r','Linewidth',line_width)
%  plot(f_axis,power_syn(:,1),'g','Linewidth',line_width)
%  plot(f_axis,power_Csyn(:,1),'b','Linewidth',line_width)
%  plot(f_axis,power_invsyn(:,1),'m','Linewidth',line_width)
% 
%   axis([0 100 -15 0])
%  legend('noise power spectra','synthetic data spectra','conditioned data spectra','inversion data spectra','Location','SouthEast')
%  title('Plots of the frequency spectra of the estimated data','FontSize',title_fontsize)
%  
%  subplot('position',[0.1 0.625 0.8 0.12]) 
%   hold on
%  plot(f_axis,power_noise(:,2),'r','Linewidth',line_width)
%  plot(f_axis,power_syn(:,2),'g','Linewidth',line_width)
%  plot(f_axis,power_Csyn(:,2),'b','Linewidth',line_width)
%  plot(f_axis,power_invsyn(:,2),'m','Linewidth',line_width)
%  axis([0 100 -15 0])
% 
%  subplot('position',[0.1 0.500 0.8 0.12])
%   hold on
%  plot(f_axis,power_noise(:,3),'r','Linewidth',line_width)
%  plot(f_axis,power_syn(:,3),'g','Linewidth',line_width)
%  plot(f_axis,power_Csyn(:,3),'b','Linewidth',line_width)
%  plot(f_axis,power_invsyn(:,3),'m','Linewidth',line_width) 
%  axis([0 100 -15 0])
%  
%  subplot('position',[0.1 0.375 0.8 0.12])
%  hold on
% 
%  plot(f_axis,power_noise(:,4),'r','Linewidth',line_width)
%  plot(f_axis,power_syn(:,4),'g','Linewidth',line_width)
%  plot(f_axis,power_Csyn(:,4),'b','Linewidth',line_width)
%  plot(f_axis,power_invsyn(:,4),'m','Linewidth',line_width) 
%  axis([0 100 -15 0])
% ylabel('dB relative to maximum value','FontSize',axis_title_fontsize)
%  subplot('position',[0.1 0.250 0.8 0.12])
%   hold on
%  plot(f_axis,power_noise(:,5),'r','Linewidth',line_width)
%  plot(f_axis,power_syn(:,5),'g','Linewidth',line_width)
%  plot(f_axis,power_Csyn(:,5),'b','Linewidth',line_width)
%  plot(f_axis,power_invsyn(:,5),'m','Linewidth',line_width)
%  axis([0 100 -15 0])
%  subplot('position',[0.1 0.125 0.8 0.12])
%  hold on
%  plot(f_axis,power_noise(:,6),'r','Linewidth',line_width)
%  plot(f_axis,power_syn(:,6),'g','Linewidth',line_width)
%  plot(f_axis,power_Csyn(:,6),'b','Linewidth',line_width)
%  plot(f_axis,power_invsyn(:,6),'m','Linewidth',line_width)
%  axis([0 100 -15 0])
%  xlabel('Frequency (Hz)','FontSize',axis_title_fontsize)
%   set(gca,'FontSize',axis_number_fontsize)

%  figure(9)
%  set(9,'Units','inches','Position',[0 0 23 16]);
%  hold on
%  title('Correlation of estimates with original data, compared to correlation with noise','FontSize',axis_title_fontsize)
%  plot(angles,corr_Csyn,'xb--','Linewidth',line_width2,'Markersize',marker_size)
%  plot(angles,corr_inv,'bo--','Linewidth',line_width2,'Markersize',marker_size)
%  plot(angles,corr_Csyn_n,'xk-','Linewidth',line_width2,'Markersize',marker_size)
%  plot(angles,corr_inv_n,'ko-','Linewidth',line_width2,'Markersize',marker_size)
%  set(gca,'FontSize',axis_number_fontsize)
%  legend('Conditioned estimate with synthetic','DIGI estimate with synthetic','Cond. with noise','DIGI with noise')
% ylabel('Angle','FontSize',axis_title_fontsize)

% do similarity? cross correlation. see whether the amplitude spectrum of
% the conditioned data is more similar to the data or the noise, see how
% much it affects

%%

% these plots show the frequency spectra of the wavelets, and will plot the
% spectrum of the noise series if coloured nosie is turned on

% figure (8)
% 
%  if coloured_noise==1
%      k=1;
%  else k=0;
%  end
% 
%  subplot(6,1,1)
%  hold on
%  plot(f_rw,power_real_w(1:length(f_rw),1),'r')
%  plot(f_c,power_w_est_C(1:length(f_c),1),'g')
%  plot(f_inv,power_w_est_inv(1:length(f_inv),1),'b')
%   axis([0 125 -25 0])
%   if k==1
%  plot(f_axis,power_noise(:,1),'m')
%  end
%  legend('real wavelet','conditioned estimated wavelet','inversion estimated wavelet','noise wavelet')
%  title('plots of the frequency spectra of the wavelets, original and estimated')
%  subplot(6,1,2) 
%   hold on
%  plot(f_rw,power_real_w(1:length(f_rw),2),'r')
%  plot(f_c,power_w_est_C(1:length(f_c),2),'g')
%  plot(f_inv,power_w_est_inv(1:length(f_inv),2),'b')
%    axis([0 125 -25 0])
%    if k==1
%  plot(f_axis,power_noise(:,2),'m')
%  end
%  subplot(6,1,3) 
%   hold on
%  plot(f_rw,power_real_w(1:length(f_rw),3),'r')
%  plot(f_c,power_w_est_C(1:length(f_c),3),'g')
%  plot(f_inv,power_w_est_inv(1:length(f_inv),3),'b')
%     axis([0 125 -25 0])
%     if k==1
%  plot(f_axis,power_noise(:,3),'m')
%  end
%  subplot(6,1,4) 
%   hold on
%   plot(f_rw,power_real_w(1:length(f_rw),4),'r')
%  plot(f_c,power_w_est_C(1:length(f_c),4),'g')
%  plot(f_inv,power_w_est_inv(1:length(f_inv),4),'b')
%    axis([0 125 -25 0])
%   if k==1
%  plot(f_axis,power_noise(:,4),'m')
%  end
%  subplot(6,1,5) 
%   hold on
%  plot(f_rw,power_real_w(1:length(f_rw),5),'r')
%  plot(f_c,power_w_est_C(1:length(f_c),5),'g')
%  plot(f_inv,power_w_est_inv(1:length(f_inv),5),'b')
%   axis([0 125 -25 0])
%  if k==1
%  plot(f_axis,power_noise(:,5),'m')
%  end
%  subplot(6,1,6) 
%   hold on
%   plot(f_rw,power_real_w(1:length(f_rw),6),'r')
%  plot(f_c,power_w_est_C(1:length(f_c),6),'g')
%   axis([0 125 -25 0])
%    plot(f_inv,power_w_est_inv(1:length(f_inv),6),'b')
%  if k==1
%  plot(f_axis,power_noise(:,6),'m')
%  end
