% Make synthetic I and G
I = zeros(1000,1);
G = zeros(1000,1);
for i = 1:1000
    if randn > 0.9
        if randn > 0.1
            I(i) = -0.1 + (0.2)*rand + 0.03*randn;
            G(i) = I(i) + 0.03*randn;
        else
            I(i) = -0.25 + (0.5)*rand + 0.03*randn;
            G(i) = -I(i)+ 0.03*randn;
        end
    end
end

% Make wavelets
%freqs = [25,20,15,10];

% for i=1:length(freqs)
%     w(:,i) = ricker_wavelet_function(freqs(i));
% end

% Estimate wavelets
Sft = abs(fft(S));
wraw = ifft(Sft/sqrt(length(I)),'symmetric');
for i=1:length(angles)
    wtmp = (2/sqrt(pi))*wraw(1:51,i).*cos(.5*pi*(0:50)'/(50.5));
    sn0 = wtmp(1)+2*sum(wtmp(2:end).*(-1).^(1:50)');
    wtmp(1) = wtmp(1)-sn0;
    wtmp = [wtmp(51:-1:2); wtmp];
    wtmp = wtmp/(norm(wtmp)/sqrt(length(wtmp)));
    w(:,i) = wtmp;
end

% Make angles
angles = [5,15,25,35];
sin2theta = (sin(angles.*pi/180)).*(sin(angles.*pi/180));

% Make reflection coeffs
for i=1:length(angles)
    R(:,i) = I+G*sin2theta(i);
end

% Make true synthetic
for i=1:length(angles)
    S(:,i) = conv2(R(:,i),w(:,i),'same');
end

% Make 'conditioned' synthetic
for i=1:length(angles)
    Sc(:,i) = conv2(R(:,i),w(:,2),'same');
end

% Estiamte I and G from not conditioned synthetic
for i=1:length(I)
    tmp = robustfit(sin2theta,S(i,:),'ols');
    Is(i,1) = tmp(1,1);
    Gs(i,1) = tmp(2,1);
end

% Estiamte I and G from conditioned synthetic
for i=1:length(I)
    tmp = robustfit(sin2theta,Sc(i,:),'ols');
    Isc(i,1) = tmp(1,1);
    Gsc(i,1) = tmp(2,1);
end

% Estimate wavelet scalers
for i=1:length(angles)
    %scale(1,i) = (norm(S)/(norm(R)*norm(w)));
    scale = [1,1,1,1];
end

% Make wavelet convolution matrix
for i = 1:length(angles)
    W(:,length(I)*(i-1)+1:length(I)*i) = convmtx(scale(i)*w(:,i),length(I));
end
clip = (size(W,1)-length(I))/2;
W = W(clip+1:end-clip,:);

% Make whitenning matrix
tmp_white = eye(length(I),length(I));
white = repmat(tmp_white,length(angles),2);
lambdaI = 0;
lambdaG = 0;
white(:,1:length(I)) = white(:,1:length(I))*lambdaI;
white(:,1+length(I):end) = white(:,1+length(I):end)*lambdaG;

% Make operator without linear trend
op = [W',W'];
for i=1:length(angles)
    op(1+(i-1)*length(I):i*length(I),1+length(I):end) = op(1+(i-1)*length(I):i*length(I),1+length(I):end)*sin2theta(i);
end
op = op+white;

% Make operater with linear trend
op_linfit = op;
c = tan(45*pi/180);
for i=1:length(angles)
    op_linfit(1+(i-1)*length(I):i*length(I),1:length(I)) = op_linfit(1+(i-1)*length(I):i*length(I),1:length(I))*(1-(sin2theta(i)/c));
end
op_linfit = op_linfit+white;

% Do inversions
inv = lsqr(op,reshape(S,[],1),1e-6,1000,[],[],[S(:,1);-S(:,1)]);
inv_linfit = lsqr(op_linfit,reshape(S,[],1),1e-6,1000,[],[],[S(:,1);zeros(length(I),1)]);

Ii = inv(1:length(I),1);
Gi = inv(1+length(I):end,1);

Iilf = inv_linfit(1:length(I),1);
Gilf = inv_linfit(1+length(I):end,1)-(Iilf/c);

% Threshold to make sparser solution
Iizscore = zscore(Ii);
Gizscore = zscore(Gi);
thresh = 0.5;
for i=1:length(I)
    if (abs(Iizscore(i)) < thresh) && (abs(Gizscore(i)) < thresh)
        Ii_clip(i) = (abs(Iizscore(i))/thresh)*Ii(i);
        Gi_clip(i) = (abs(Gizscore(i))/thresh)*Gi(i);
    else
        Ii_clip(i) = Ii(i);
        Gi_clip(i) = Gi(i);
    end
end

Iilfzscore = zscore(Iilf);
Gilfzscore = zscore(Gilf);
thresh = 0.5;
for i=1:length(I)
    if (abs(Iilfzscore(i)) < thresh) && (abs(Gilfzscore(i)) < thresh)
        Iilf_clip(i) = (abs(Iilfzscore(i))/thresh)*Iilf(i);
        Gilf_clip(i) = (abs(Gilfzscore(i))/thresh)*Gilf(i);
    else
        Iilf_clip(i) = Iilf(i);
        Gilf_clip(i) = Gilf(i);
    end
end


