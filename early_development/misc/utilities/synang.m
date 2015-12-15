clear all

%% make synthetic

nt = 1000;
sparsity = 0.7; % (0-1) 1 = fully sparse
max_amplitude = 0.1;
ricker_cf = [37 34 31 28 25 22 19 16];
ricker_sampling = 0.001;
background_dominace = 0.8; % (0-1) 1 = all background
chi = 18;
chinoise = 2;
IGnoise = 0.8; % (0-1) 0 = no noise
angles = [5 10 15 20 25 30 35 40];
subsample = [1;zeros((0.004/ricker_sampling)-1,1)];
subsample = bsxfun(@times,ones(1,floor(nt/length(subsample))),subsample);
subsample = logical(subsample(:));
anomaly_multiplier = 2; % (>0)

ref = synref(nt,sparsity,max_amplitude);
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

Gnoise = IGnoise*std(G)*randn(1000,1);
G(~isnan(class)) = G(~isnan(class))+Gnoise(~isnan(class));

Inoise = IGnoise*std(I)*randn(1000,1);
I(~isnan(class)) = I(~isnan(class))+Inoise(~isnan(class));

eref = bsxfun(@plus,I,bsxfun(@times,G,sin(angles*pi/180).*sin(angles*pi/180)));

for ii = 1:length(angles)
    w_tmp{ii} = ricker(ricker_cf(ii),ricker_sampling);
end

max_w = length(w_tmp{length(angles)});

for ii = 1:length(angles)
    w_tmp{ii} = [zeros((max_w-length(w_tmp{ii}))/2,1);w_tmp{ii};zeros((max_w-length(w_tmp{ii}))/2,1)];
    w(:,ii) = w_tmp{ii};
end

for ii = 1:length(angles)
    syn(:,ii) = conv2(eref(:,ii),w(:,ii),'same');
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

% West = abs(fft(syn));
% West = West(1:nt/2,:);
% West = bsxfun(@times,West,(cos(0:pi/(2*((nt/2)-1)):pi/2))');
% West = circshift(ifft(West,'symmetric'),nt/4);
% West = bsxfun(@times,West,[cos(-pi/2:pi/(2*((nt/4)-1)):0)';ones(250,1)]);
% West = bsxfun(@times,West,[ones(250,1);cos(0:pi/(2*((nt/4)-1)):pi/2)']);

for ii = 1:length(angles)
    filter(:,ii) = match_call(w(:,length(angles)),w(:,ii),-6);
    %filter(:,ii) = match_call(West(:,4),West(:,ii),-6);
end

for ii = 1:length(angles)
    Csyn(:,ii) = conv2(syn(:,ii),filter(:,ii),'same');
end

Cava = [ones(size(angles')) sin(angles'*pi/180).*sin(angles'*pi/180)]\Csyn';

CI = Cava(1,:)';
CG = Cava(2,:)';
CIG = CI.*CG;
CIG(CIG<0) = 0;

%% Inverted AVA
Iwhite = 0;
Gwhite = 0;
alpha = 1;
wsmooth = 10;

for ii = 1:length(angles)
    Iw(:,ii) = w(:,ii) + [zeros(length(w)/2,1);Iwhite;zeros(length(w)/2,1)];
    Gw(:,ii) = w(:,ii)*(sin(angles(ii)*pi/180)*sin(angles(ii)*pi/180)) + [zeros(length(w)/2,1);Gwhite;zeros(length(w)/2,1)];
    
    Imatrix_tmp = sparse(convmtx(Iw(:,ii),nt));
    Imatrix_tmp = Imatrix_tmp(length(w)/2:end-length(w)/2,:);
    
    Gmatrix_tmp = sparse(convmtx(Gw(:,ii),nt));
    Gmatrix_tmp = Gmatrix_tmp(length(w)/2:end-length(w)/2,:);
    
    if ii == 1
        Imatrix = Imatrix_tmp;
        Gmatrix = Gmatrix_tmp;
    else
        Imatrix = sparse([Imatrix; Imatrix_tmp]);
        Gmatrix = sparse([Gmatrix; Gmatrix_tmp]);
    end
end

for ii = 1:length(angles)
    synI(:,ii) = conv2(CI,w(:,ii),'same');
    synIw(:,ii) = conv2(CI*sin(angles(ii)*pi/180)*sin(angles(ii)*pi/180),w(:,ii),'same');
end

chiest = 180/pi*atan(norm(synIw)/norm(synI-syn));

EERmatrix = sparse([cos(chiest*pi/180)*convmtx(w(:,7),nt) sin(chiest*pi/180)*convmtx(w(:,7),nt)]);
EERmatrix = alpha*EERmatrix(length(w)/2:end-length(w)/2,:);

smooth = spdiags([-wsmooth*ones(2*nt,1) 2*wsmooth*ones(2*nt,1) -wsmooth*ones(2*nt,1)],[-1 0 1],2*nt,2*nt);

IGmatrix = sparse([Imatrix Gmatrix; EERmatrix; smooth]);

% [spikes ~] = sparse_decon(syn(:,7),w(:,7),0.5,500,0);

data = [syn(:);syn(:,7);zeros(2*nt,1)];
% model = [syn(:,1);-syn(:,1)/tan(chi*180/pi)];
model = [CI;CG];

% itermax = 20;
% sc = 0.0001;
% mu = 0.0001;
% squraeIGmatrix = sparse(IGmatrix'*IGmatrix);
% R0 = trace(squraeIGmatrix);
% Q = sparse(mu*R0*eye(2*nt));
% 
% for kk = 1:itermax
%     g = IGmatrix'*data;
%     wIGmatrix = squraeIGmatrix + Q;
%     invava = wIGmatrix\g;
%     Q = sparse(mu*diag(1./(abs(invava)+sc)));
% end

invava = lsqr(IGmatrix,data,1e-3,500,[],[],model);

invI = invava(1:nt);
invG = invava(nt+1:end);

% index = or(abs(invI)<5e-3,abs(invG)<5e-3);
% invI(index) = 0;
% invG(index) = 0;

% invIzscore = zscore(invI);
% invGzscore = zscore(invG);
% thresh = 0;
% for i=1:length(invI)
%     if (abs(invIzscore(i)) < thresh) && (abs(invGzscore(i)) < thresh)
%         invI(i) = 0;%(abs(invIzscore(i))/thresh)*invI(i);
%         invG(i) = 0;%(abs(invGzscore(i))/thresh)*invG(i);
%     else
%         invI(i) = invI(i);
%         invG(i) = invG(i);
%     end
% end

invIG = invI.*invG;
invIG(invIG<0) = 0;


%% plots

scaledNCI = NCI*norm(I)/norm(NCI(subsample));
scaledNCG = NCG*norm(G)/norm(NCG(subsample));
scaledCI = CI*norm(I)/norm(CI(subsample));
scaledCG = CG*norm(G)/norm(CG(subsample));
scaledinvI = invI*norm(I)/norm(invI(subsample));
scaledinvG = invG*norm(G)/norm(invG(subsample));

figure(1)

subplot(3,1,1)
plot(syn)
axis tight
title('Input Synthetic Angle Traces')
xlabel('Sample number')
ylabel('Amplitude')

subplot(3,1,2)
stem(I,'marker','none')
axis tight
title('Input Synthetic Intercept')
xlabel('Sample number')
ylabel('Amplitude')

subplot(3,1,3)
stem(G,'marker','none')
axis tight
title('Input Synthetic Gradient')
xlabel('Sample number')
ylabel('Amplitude')


figure(2)

subplot(3,2,1)
stem(I,'marker','none','color',[150/255 150/255 150/255])
hold all
plot(scaledNCI)
hold off
axis tight
title('Intercept from raw amplitude regression (truth in grey)')
xlabel('Sample number')
ylabel('Amplitude')

subplot(3,2,2)
stem(G,'marker','none','color',[150/255 150/255 150/255])
hold all
plot(scaledNCG)
hold off
axis tight
title('Gradient from raw amplitude regression (truth in grey)')
xlabel('Sample number')
ylabel('Amplitude')

subplot(3,2,3)
stem(I,'marker','none','color',[150/255 150/255 150/255])
hold all
plot(scaledCI)
hold off
axis tight
title('Intercept from conditioned amplitude regression (truth in grey)')
xlabel('Sample number')
ylabel('Amplitude')

subplot(3,2,4)
stem(G,'marker','none','color',[150/255 150/255 150/255])
hold all
plot(scaledCG)
hold off
axis tight
title('Gradient from conditioned amplitude regression (truth in grey)')
xlabel('Sample number')
ylabel('Amplitude')

subplot(3,2,5)
stem(I,'marker','none','color',[150/255 150/255 150/255])
hold all
plot(scaledinvI)
hold off
axis tight
title('Intercept from inversion (truth in grey)')
xlabel('Sample number')
ylabel('Amplitude')

subplot(3,2,6)
stem(G,'marker','none','color',[150/255 150/255 150/255])
hold all
plot(scaledinvG)
hold off
axis tight
title('Gradient from inversion (truth in grey)')
xlabel('Sample number')
ylabel('Amplitude')


figure(3)

bpI = conv(I,w(:,1),'same');
bpG = conv(G,w(:,1),'same');

subplot(2,2,1)
time = (1:1:nt);
% time(IG==0) = 0;
scatter(bpI(subsample),bpG(subsample),[],time(subsample),'filled');
axis([-0.5 0.5 -0.5 0.5]);
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Sample Number');
title('True Intercept vs Gradient convolved filtered to seismic bandwidth')
xlabel('Intercept')
ylabel('Gradient')

subplot(2,2,2)
% time = (1:1:nt);
% time(NCIG==0) = 0;
% scatter(NCI(subsample),NCG(subsample),[],time(subsample),'filled'); axis equal
colours = ((scaledNCI - conv(I,w(:,1),'same')).^2 + (scaledNCG - conv(G,w(:,1),'same')).^2).^(0.5);
scatter(scaledNCI(subsample),scaledNCG(subsample),[],colours(subsample),'filled');
axis([-0.5 0.5 -0.5 0.5]);
caxis([0 0.3]);
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Euclidean distance from true IG location');
title('Intercept vs Gradient from raw amplitude regression')
xlabel('Intercept')
ylabel('Gradient')

subplot(2,2,3)
% time = (1:1:nt);
% time(CIG==0) = 0;
% scatter(CI(subsample),CG(subsample),[],time(subsample),'filled'); axis equal
colours = ((scaledCI - conv(I,w(:,1),'same')).^2 + (scaledCG - conv(G,w(:,1),'same')).^2).^(0.5);
scatter(scaledCI(subsample),scaledCG(subsample),[],colours(subsample),'filled');
axis([-0.5 0.5 -0.5 0.5]);
caxis([0 0.3]);
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Euclidean distance from true IG location');
title('Intercept vs Gradient from conditioned amplitude regression')
xlabel('Intercept')
ylabel('Gradient')

subplot(2,2,4)
% time = (1:1:nt);
% time(invIG==0) = 0;
% scatter(invI(subsample),invG(subsample),[],time(subsample),'filled'); axis equal
colours = ((scaledinvI - conv(I,w(:,1),'same')).^2 + (scaledinvG - conv(G,w(:,1),'same')).^2).^(0.5);
scatter(scaledinvI(subsample),scaledinvG(subsample),[],colours(subsample),'filled');
axis([-0.5 0.5 -0.5 0.5]);
caxis([0 0.3]);
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Euclidean distance from true IG location');
title('Intercept vs Gradient from inversion')
xlabel('Intercept')
ylabel('Gradient')


figure(4)

subplot(3,1,1)
stem(I*cos(chi*pi/180) + G*sin(chi*pi/180),'marker','none','color',[150/255 150/255 150/255])
hold all
plot(scaledNCI*cos(chi*pi/180) + scaledNCG*sin(chi*pi/180))
hold off
axis tight
title('Minimum energy EER projection using I and G from raw amplitude regression (truth in grey)')
xlabel('Sample Number')
ylabel('Amplitude')

subplot(3,1,2)
stem(I*cos(chi*pi/180) + G*sin(chi*pi/180),'marker','none','color',[150/255 150/255 150/255])
hold all
plot(scaledCI*cos(chi*pi/180) + scaledCG*sin(chi*pi/180))
hold off
axis tight
title('Minimum energy EER projection using I and G from conditioned amplitude regression (truth in grey)')
xlabel('Sample Number')
ylabel('Amplitude')

subplot(3,1,3)
stem(I*cos(chi*pi/180) + G*sin(chi*pi/180),'marker','none','color',[150/255 150/255 150/255])
hold all
plot(scaledinvI*cos(chi*pi/180) + scaledinvG*sin(chi*pi/180))
hold off
axis tight
title('Minimum energy EER projection using I and G from inversion (truth in grey)')
xlabel('Sample Number')
ylabel('Amplitude')


figure(5)

ncolours = 64;
colours = [[ones(ncolours/2,1);(1:-1/((ncolours/2)-1):0)'] [(0:1/((ncolours/2)-2):1)';1;1;(1:-1/((ncolours/2)-2):0)'] [(0:1/((ncolours/2)-1):1)';ones(ncolours/2,1)]];

for ii = 1:length(angles)
    NCsyn(:,ii) = conv(NCI + NCG*(sin(angles(ii)*pi/180)*sin(angles(ii)*pi/180)),w(:,ii),'same');
    Csyn(:,ii) = conv(CI + CG*(sin(angles(ii)*pi/180)*sin(angles(ii)*pi/180)),w(:,ii),'same');
    invsyn(:,ii) = conv(invI + invG*(sin(angles(ii)*pi/180)*sin(angles(ii)*pi/180)),w(:,ii),'same');
end

subplot(3,3,1)
imagesc(angles,(1:1:nt),syn,[-0.9*max(max(abs(syn))) 0.9*max(max(abs(syn)))])
axis fill
colormap(colours);
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Amplitude');
title('(1) Input Synthetic Angle Traces')
xlabel('Angles (degrees)')
ylabel('Sample Number')

subplot(3,3,2)
imagesc(angles,(1:1:nt),NCsyn*norm(syn)/norm(NCsyn),[-0.9*max(max(abs(syn))) 0.9*max(max(abs(syn)))])
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Amplitude');
title({'(2) Angle Traces forward modelled using';'I and G from raw amplitude regression'})
xlabel('Angles (degrees)')
ylabel('Sample Number')

subplot(3,3,3)
imagesc(angles,(1:1:nt),syn - (NCsyn*norm(syn)/norm(NCsyn)),[-0.9*max(max(abs(syn))) 0.9*max(max(abs(syn)))])
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Amplitude');
title('Residual of (1) and (2)')
xlabel('Angles (degrees)')
ylabel('Sample Number')

subplot(3,3,4)
imagesc(angles,(1:1:nt),syn,[-0.9*max(max(abs(syn))) 0.9*max(max(abs(syn)))])
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Amplitude');
title('(1) Input Synthetic Angle Traces')
xlabel('Angles (degrees)')
ylabel('Sample Number')

subplot(3,3,5)
imagesc(angles,(1:1:nt),Csyn*norm(syn)/norm(Csyn),[-0.9*max(max(abs(syn))) 0.9*max(max(abs(syn)))])
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Amplitude');
title({'(3) Angle Traces forward modelled using';'I and G from conditioned amplitude regression'})
xlabel('Angles (degrees)')
ylabel('Sample Number')

subplot(3,3,6)
imagesc(angles,(1:1:nt),syn - (Csyn*norm(syn)/norm(Csyn)), [-0.9*max(max(abs(syn))) 0.9*max(max(abs(syn)))])
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Amplitude');
title('Residual of (1) and (3)')
xlabel('Angles (degrees)')
ylabel('Sample Number')

subplot(3,3,7)
imagesc(angles,(1:1:nt),syn,[-0.9*max(max(abs(syn))) 0.9*max(max(abs(syn)))])
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Amplitude');
title('(1) Input Synthetic Angle Traces')
xlabel('Angles (degrees)')
ylabel('Sample Number')

subplot(3,3,8)
imagesc(angles,(1:1:nt),invsyn,[-0.9*max(max(abs(invsyn))) 0.9*max(max(abs(invsyn)))])
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Amplitude');
title('(4) Angle Traces forward modelled using I and G from inversion')
xlabel('Angles (degrees)')
ylabel('Sample Number')

subplot(3,3,9)
imagesc(angles,(1:1:nt),syn - (invsyn*norm(syn)/norm(invsyn)), [-0.9*max(max(abs(syn))) 0.9*max(max(abs(syn)))])
colorbar
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Amplitude');
title('Residual of (1) and (4)')
xlabel('Angles (degrees)')
ylabel('Sample Number')



%%
% %% non-linear inversion
% % model = [CI;CG;2.7*ones(nt,1)];
% trend = (1/3)*ones(nt,1);
% trend(CIG>0) = -3;
% model = [invI;invG;trend];
% itermax = 100;
% % data = [syn(:);zeros(nt,1);2.7*ones(nt,1)];
% data = [syn(:);zeros(nt,1);trend;zeros(4*nt,1)];
% Ismooth = 100;
% Gsmooth = 100;
% Idamp = 10;
% Gdamp = 10;
% 
% nlIGmatrix = [Imatrix Gmatrix zeros(nt*length(angles),nt); eye(nt) zeros(nt,nt) diag(model(nt+1:2*nt,1)); zeros(nt,nt) zeros(nt,nt) eye(nt)];
% jacobian = [Imatrix Gmatrix zeros(nt*length(angles),nt); eye(nt) diag(model(2*nt+1:end,1)) diag(model(nt+1:2*nt,1)); zeros(nt,nt) zeros(nt,nt) eye(nt)];
% damp = [diag([(Idamp)*ones(nt,1); (Gdamp)*ones(nt,1)]) zeros(2*nt,nt)];
% smooth = [spdiags([[(-Ismooth/2)*ones(nt,1); (-Gsmooth/2)*ones(nt,1)] [Ismooth*ones(nt,1); Gsmooth*ones(nt,1)] [(-Ismooth/2)*ones(nt,1); (-Gsmooth/2)*ones(nt,1)]],[-1 0 1],2*nt,2*nt) zeros(2*nt,nt)];
% 
% nlIGmatrix = [nlIGmatrix; damp; smooth];
% jacobian = [jacobian; damp; smooth];
% 
% residuals = (nlIGmatrix*model - data);
% 
% for kk = 1:20
%     
%     delta = (jacobian'*jacobian)\(jacobian'*residuals);
% %     delta = lsqr(jacobian,residuals,1e-3,500);
%     
%     if kk == 1
%         invava = model + delta;
%     else
%         invava = invava + delta;
%     end
%     
%     residuals = (nlIGmatrix*invava - data);
%     jacobian = [Imatrix Gmatrix zeros(nt*length(angles),nt); eye(nt) diag(invava(2*nt+1:end,1)) diag(invava(nt+1:2*nt,1)); zeros(nt,nt) zeros(nt,nt) eye(nt)];
%     jacobian = [jacobian; damp; smooth];
% end
%     
% % 
% % % objective = (nlIGmatrix*model - data).^2;
% % % jacopModel = jacobian'*objective;
% % % change = ((jacobian*model) - data).^2;
% % 
% % % g_old = 2*nlIGmatrix'*((nlIGmatrix*model)-data);
% % % g_init = -(jacobian'*jacobian*model - jacobian'*data);
% % % g_init = -(jacobian'*jacobian*model);
% % % g_init = -(2*nlIGmatrix'*((nlIGmatrix*model)-data));
% % g_init = -(2*nlIGmatrix'*((nlIGmatrix*model)-data));
% % 
% % step = 0.1;
% % res_step(1,1) = norm(nlIGmatrix*model - data);
% % res_step(1,2) = norm(nlIGmatrix*(model + step.*g_init) - data);
% % while res_step(1,2) > res_step(1,1)
% %     step = 0.1*step;
% %     res_step(1,2) = norm(nlIGmatrix*(model + step.*g_init) - data);
% %     if step < 1e-12
% %         break
% %     end
% % end
% % 
% % % step = 0.1*ones(3*nt,1);
% % % res_step(:,1) = (nlIGmatrix'*nlIGmatrix*model - nlIGmatrix'*data).^2;
% % % res_step(:,2) = (nlIGmatrix'*nlIGmatrix*(model + step.*g_init) - nlIGmatrix'*data).^2;
% % % sqnlIGmatrix = nlIGmatrix'*nlIGmatrix;
% % % for ii = 1:length(step)
% % %     while res_step(ii,2) > res_step(ii,1)
% % %         step(ii,1) = 0.1*step(ii,1);
% % %         res_step(ii,2) = (sqnlIGmatrix(ii,:)*(model + step(ii,1)*g_init) - nlIGmatrix(:,ii)'*data).^2;
% % %         if step(ii,1) < 1e-6
% % %             step(ii,1) = 0;
% % %             break
% % %         end
% % %     end
% % % end
% % 
% % invava = model + step.*g_init;
% % 
% % stop = 0;
% % 
% % for kk = 1:itermax
% %     residual(kk,1) = norm(nlIGmatrix*(invava) - data);    
% % %     jacobian = [Imatrix Gmatrix zeros(nt*length(angles),nt); eye(nt) diag(invava(2*nt+1:end,1)) diag(invava(nt+1:2*nt,1)); zeros(nt,nt) zeros(nt,nt) eye(nt)];
% %     nlIGmatrix = [Imatrix Gmatrix zeros(nt*length(angles),nt); eye(nt) zeros(nt,nt) diag(invava(nt+1:2*nt,1)); zeros(nt,nt) zeros(nt,nt) eye(nt)];
% %     if kk == 1
% % %         g = -(jacobian'*jacobian*invava - jacobian'*data);
% % %         g = -(jacobian'*jacobian*invava);
% %         g = -(2*nlIGmatrix'*((nlIGmatrix*model)-data));
% %         B = (g'*g)/(g_init'*g_init);
% %         S = g + B*g_init;
% %         step = 0.1;
% %         res_step(kk,1) = norm(nlIGmatrix*invava - data);
% %         res_step(kk,2) = norm(nlIGmatrix*(invava + step.*S) - data);
% %         while res_step(kk,2) > res_step(kk,1)
% %             step = 0.1*step;
% %             res_step(1,2) = norm(nlIGmatrix*(invava + step.*S) - data);
% %             if step < 1e-12
% %                 stop = 1;
% %                 break
% %             end
% %         end
% %         if stop == 1
% %             break
% %         end
% %         invava = invava + step*S;
% %     else
% %         g_old = g;
% % %         g = -(jacobian'*jacobian*invava - jacobian'*data);
% % %         g = -(jacobian'*jacobian*invava);
% %         g = -(2*nlIGmatrix'*((nlIGmatrix*model)-data));
% %         B = (g'*g)/(g_old'*g_old);
% %         S = g + B*S;
% %         step = 0.1;
% %         res_step(kk,1) = norm(nlIGmatrix*invava - data);
% %         res_step(kk,2) = norm(nlIGmatrix*(invava + step.*S) - data);
% %         while res_step(kk,2) > res_step(kk,1)
% %             step = 0.1*step;
% %             res_step(1,2) = norm(nlIGmatrix*(invava + step.*S) - data);
% %             if step < 1e-12
% %                 stop = 1;
% %                 break
% %             end
% %         end
% %         if stop == 1
% %             break
% %         end
% %         invava = invava + step*S;
% %     end
% % end
% % 
% % % invava = model - step*(jacopModel);
% % 
% % % for kk = % 1:itermax
% % %     residual(kk,1) = (sum(objective))/length(objective);
% % %     if residual < 1e-4
% % %         break
% % %     end
% % %     
% % %     nlIGmatrix = [Imatrix Gmatrix zeros(nt*length(angles),nt); eye(nt) zeros(nt,nt) diag(invava(nt+1:2*nt,1)); zeros(nt,nt) zeros(nt,nt) eye(nt)];
% % %     objective = (nlIGmatrix*invava - data).^2;
% % % 
% % %     jacobian = [Imatrix Gmatrix zeros(nt*length(angles),nt); eye(nt) diag(invava(2*nt+1:end,1)) diag(invava(nt+1:2*nt,1)); zeros(nt,nt) zeros(nt,nt) eye(nt)];
% % %     jacopModel = jacobian'*objective;
% % % 
% % %     step = 0.0001;
% % % 
% % %     invava = model - step*(jacopModel);
% % % end
% %   
% magicI = invava(1:nt);
% magicG = invava(nt+1:2*nt);
% cosmo = invava(2*nt+1:end);
% 
% magicIG = magicI.*magicG;
% magicIG(magicIG<0) = 0;

%%

% nlChi = 18;
% data = syn(:);
% model = [CI;zeros(nt,1);zeros(nt,1);nlChi];
% 
% for kk = 1:5
% 
%     if kk == 1
%         nlI = model(1:nt,1);
%         nldI = model(nt+1:2*nt,1);
%         nldG = model(2*nt+1:3*nt,1);
%     end
% 
%     cosChi = cos(nlChi*pi/180);
%     sinChi = sin(nlChi*pi/180);
%     tanChi = tan(nlChi*pi/180);
% 
%     for ii = 1:length(angles)
%         blockI(:,ii) = w(:,ii)*(1-((sin(angles(ii)*pi/180)*sin(angles(ii)*pi/180))/(tanChi)));
%         blockdI(:,ii) = w(:,ii)*(cosChi);
%         blockdG(:,ii) = w(:,ii)*(sinChi);
%         blockChi(:,ii) = (sin(angles(ii)*pi/180)*sin(angles(ii)*pi/180))/(tanChi^2*cosChi^2) - nldI*sinChi + nldG.*cosChi;
% 
%         blockopdI(:,ii) = w(:,ii)*(sinChi);
%         blockopdG(:,ii) = w(:,ii)*(cosChi);
% 
%         blockI_tmp = convmtx(blockI(:,ii),nt);
%         blockI_tmp = blockI_tmp(length(w)/2:end-length(w)/2,:);
% 
%         blockdI_tmp = convmtx(blockdI(:,ii),nt);
%         blockdI_tmp = blockdI_tmp(length(w)/2:end-length(w)/2,:);
% 
%         blockdG_tmp = convmtx(blockdI(:,ii),nt);
%         blockdG_tmp = blockdG_tmp(length(w)/2:end-length(w)/2,:);
% 
%         blockChi_tmp = blockChi(:,ii);
% 
%         blockopdI_tmp = convmtx(blockopdI(:,ii),nt);
%         blockopdI_tmp = blockopdI_tmp(length(w)/2:end-length(w)/2,:);
% 
%         blockopdG_tmp = convmtx(blockopdG(:,ii),nt);
%         blockopdG_tmp = blockopdG_tmp(length(w)/2:end-length(w)/2,:);
% 
%         if ii == 1
%             jacI = blockI_tmp;
%             jacdI = blockdI_tmp;
%             jacdG = blockdG_tmp;
%             jacChi = blockChi_tmp;
%             opdI = blockopdI_tmp;
%             opdG = blockopdG_tmp;
%         else
%             jacI = [jacI; blockI_tmp];
%             jacdI = [jacdI; blockdI_tmp];
%             jacdG = [jacdG; blockdG_tmp];
%             jacChi = [jacChi; blockChi_tmp];
%             opdI = [opdI; blockopdI_tmp];
%             opdG = [opdG; blockopdG_tmp];
%         end
% 
%     end
% 
%     operator = sparse([jacI opdI opdG zeros(nt*length(angles),1)]);
%     jacobian = sparse([jacI jacdI jacdG jacChi]);
% 
%     if kk == 1
%         residuals = (operator*model -data);
%     else
%         residuals = (operator*invava -data);
%     end
% 
%     delta = (jacobian'*jacobian)\(jacobian'*residuals);
% 
%     if kk == 1
%         invava = model + delta;
%     else
%         invava = invava + delta;
%     end
%     nlI = invava(1:nt,1);
%     nldI = invava(nt+1:2*nt,1);
%     nldG = invava(2*nt+1:3*nt,1);
%     nlChi = invava(end,1);
%     
%     endI(:,kk) = nlI + nldI*cos(nlChi*pi/180);
%     endG(:,kk) = -nlI/tan(nlChi*pi/180) + nldG*sin(nlChi*pi/180);
% 
% end
