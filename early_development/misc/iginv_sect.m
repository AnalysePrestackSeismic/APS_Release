% Some setup parameters
I = zeros(251,1);
G = zeros(251,1);
Is = zeros(size(near));
Gs = zeros(size(near));
Ii = zeros(size(near));
Gi = zeros(size(near));
Iilf = zeros(size(near));
Gilf = zeros(size(near));
Ii_clip = zeros(size(near));
Gi_clip = zeros(size(near));
Iilf_clip = zeros(size(near));
Gilf_clip = zeros(size(near));
taper = 10; % Defines cos-taper length in samples for beginning and end of trace (might not be needed)
wrms = 1; % Scales the RMS of the wavelet to be equal to wrms
stepout = 3; % Trace stepout for inversion (total traces used equals 1+2*stepout)
residual = nan(length(near)-2*stepout);
iter = nan(length(near)-2*stepout);
chi = nan(length(near)-2*stepout);

% Estimate wavelets for near, mid and far
% nearft = mean(abs(fft(near)),2);
% wraw = ifft(nearft/sqrt(length(I)),'symmetric');
% wtmp = (2/sqrt(pi))*wraw(1:51).*cos(.5*pi*(0:50)'/(50));
% sn0 = wtmp(1)+2*sum(wtmp(2:end).*(-1).^(1:50)');
% wtmp(1) = wtmp(1)-sn0;
% wtmp = [wtmp(51:-1:2); wtmp];
% wtmp = wtmp/(norm(wtmp/wrms)/sqrt(length(wtmp)));
% wnear = wtmp;
% 
% midft = mean(abs(fft(mid)),2);
% wraw = ifft(midft/sqrt(length(I)),'symmetric');
% wtmp = (2/sqrt(pi))*wraw(1:51).*cos(.5*pi*(0:50)'/(50));
% sn0 = wtmp(1)+2*sum(wtmp(2:end).*(-1).^(1:50)');
% wtmp(1) = wtmp(1)-sn0;
% wtmp = [wtmp(51:-1:2); wtmp];
% wtmp = wtmp/(norm(wtmp/wrms)/sqrt(length(wtmp)));
% wmid = wtmp;
% 
% farft = mean(abs(fft(far)),2);
% wraw = ifft(farft/sqrt(length(I)),'symmetric');
% wtmp = (2/sqrt(pi))*wraw(1:51).*cos(.5*pi*(0:50)'/(50));
% sn0 = wtmp(1)+2*sum(wtmp(2:end).*(-1).^(1:50)');
% wtmp(1) = wtmp(1)-sn0;
% wtmp = [wtmp(51:-1:2); wtmp];
% wtmp = wtmp/(norm(wtmp/wrms)/sqrt(length(wtmp)));
% wfar = wtmp;

% w = [wnear,wmid,wfar];

% Start of inversion loop
for k = 1+stepout:length(near)-stepout
%for k = 12:12
    fprintf('\nInverting traces %d to %d of %d\n',k-stepout,k+stepout,length(near))
    
    % Taper start and end of traces
    near_taper = near(:,k-stepout:k+stepout);
    mid_taper = mid(:,k-stepout:k+stepout);
    far_taper = far(:,k-stepout:k+stepout);
    
    for i = 1:(2*stepout)+1
        near_taper(end-taper:end,i) = (2/sqrt(pi))*near(end-taper:end,k-stepout+i-1).*cos(.5*pi*(0:taper)'/taper);
        near_taper(1:taper+1,i) = (2/sqrt(pi))*near(1:taper+1,k-stepout+i-1).*cos(.5*pi*(taper:-1:0)'/taper);
        mid_taper(end-taper:end,i) = (2/sqrt(pi))*mid(end-taper:end,k-stepout+i-1).*cos(.5*pi*(0:taper)'/taper);
        mid_taper(1:taper+1,i) = (2/sqrt(pi))*mid(1:taper+1,k-stepout+i-1).*cos(.5*pi*(taper:-1:0)'/taper);
        far_taper(end-taper:end,i) = (2/sqrt(pi))*far(end-taper:end,k-stepout+i-1).*cos(.5*pi*(0:taper)'/taper);
        far_taper(1:taper+1,i) = (2/sqrt(pi))*far(1:taper+1,k-stepout+i-1).*cos(.5*pi*(taper:-1:0)'/taper);
    end
    
    S = [reshape(near_taper,[],1),reshape(mid_taper,[],1),reshape(far_taper,[],1)];
    
    % Make angles
    angles = [8,23,38];
    sin2theta = (sin(angles.*pi/180)).*(sin(angles.*pi/180));

    % Estimate wavelets
    Sft = abs(fft(S));
    wraw = ifft(Sft/sqrt(length(I)),'symmetric');
    for i=1:length(angles)
        wtmp = (2/sqrt(pi))*wraw(1:51,i).*cos(.5*pi*(0:50)'/(50.5));
        sn0 = wtmp(1)+2*sum(wtmp(2:end).*(-1).^(1:50)');
        wtmp(1) = wtmp(1)-sn0;
        wtmp = [wtmp(51:-1:2); wtmp];
        wtmp = wtmp/(norm(wtmp/wrms)/sqrt(length(wtmp)));
        w(:,i) = wtmp;
    end

    % Make reflection coeffs
    % for i=1:length(angles)
    %     R(:,i) = I+G*sin2theta(i);
    % end

    % Make true synthetic
    % for i=1:length(angles)
    %     S(:,i) = conv2(R(:,i),w(:,i),'same');
    % end

    % Make 'conditioned' synthetic
    % for i=1:length(angles)
    %     Sc(:,i) = conv2(R(:,i),w(:,2),'same');
    % end

    % Estiamte I and G from line fitting to the data
    for i=(stepout*length(I))+1:(stepout+1)*length(I)
        tmp = robustfit(sin2theta,S(i,:),'ols');
        Is(i-stepout*length(I),k) = tmp(1,1);
        Gs(i-stepout*length(I),k) = tmp(2,1);
    end

    % Estimate wavelet scalers
    % for i=1:length(angles)
        %scale(1,i) = (norm(S)/(norm(R)*norm(w)));
        scale = [1,1,1];
    % end

    % Make wavelet convolution matrix
    for i = 1:length(angles)
        W(:,length(S)*(i-1)+1:length(S)*i) = convmtx(scale(i)*w(:,i),length(S));
    end
    clip = (size(W,1)-length(S))/2;
    W = W(clip+1:end-clip,:);
    W = sparse(W);

    % Make whitenning matrix
    tmp_white = eye(length(S),length(S));
    white = repmat(tmp_white,length(angles),2);
    lambdaI = 1;
    lambdaG = 1;
    white(:,1:length(S)) = white(:,1:length(S))*lambdaI;
    white(:,1+length(S):end) = white(:,1+length(S):end)*lambdaG;
    white = sparse(white);

    % Make operator without linear trend
    op = [W',W'];
    for i=1:length(angles)
        op(1+(i-1)*length(S):i*length(S),1+length(S):end) = op(1+(i-1)*length(S):i*length(S),1+length(S):end)*sin2theta(i);
    end
    op = op+white;

    % Estimate c
    minchi = 10;
    maxchi = 60;
    carray = (tan((minchi:5:maxchi)*pi/180))';
    relres = zeros(1+(maxchi/5-minchi/5),1);
    for n=1:(maxchi/5-minchi/5)+1
        op_linfit = op;
        for i=1:length(angles)
            op_linfit(1+(i-1)*length(S):i*length(S),1:length(S)) = op_linfit(1+(i-1)*length(S):i*length(S),1:length(S))*(1-(sin2theta(i)/carray(n)));
        end
        %op_linfit = op_linfit+white;

        % Do inversion to find best c
        [~,~,relres(n,1)] = lsqr(op_linfit,reshape(S,[],1),1e-3,5,[],[],[S(:,1);zeros(length(S),1)]);
        %relrestrack(n,1) = relres;
        %relRMSdG(n,1) = (norm(inv_linfit(1+length(S):end))/sqrt(length(S)))/(norm(S)/sqrt(length(S)*length(angles)));
    end
    relresi = interpft(relres,51);
    [~,chi(k,1)] = min(relresi);
    fprintf('chi = %d\n',chi(k)+minchi-1);
    c = tan((chi(k)+minchi-1)*pi/180);
    
    % Make operater with linear trend
    op_linfit = op;
    for i=1:length(angles)
        op_linfit(1+(i-1)*length(S):i*length(S),1:length(S)) = op_linfit(1+(i-1)*length(S):i*length(S),1:length(S))*(1-(sin2theta(i)/c));
    end
    
    % Do inversions
    % inv = lsqr(op,reshape(S,[],1),1e-6,500,[],[],[S(:,1);-S(:,1)/c]);
    [inv_linfit,~,residual(k,1),iter(k,1)] = lsqr(op_linfit,reshape(S,[],1),1e-3,1000,[],[],[S(:,1);zeros(length(S),1)]);
    fprintf('convergence after %d iterations with relative residual of %f\n',iter(k),residual(k));

%     Ii(:,k) = inv(1:length(I),1);
%     Gi(:,k) = inv(1+length(I):end,1);

    % Select the kth trace as the result for the kth position
%     Iilf(:,k) = inv_linfit((stepout*length(I))+1:(stepout+1)*length(I),1);
%     Gilf(:,k) = inv_linfit(end+1-((stepout+1)*length(I)):end-(stepout*length(I)))-(Iilf(:,k)/c);

    % Average the kth +/- stepout traces for the results at the kth position
    Iilf(:,k) = mean(reshape(inv_linfit(1:length(S),1),[],1+2*stepout),2);
    Gilf(:,k) = mean(reshape(inv_linfit(1+length(S):end,1),[],1+2*stepout),2)-(Iilf(:,k)/c);

    % Threshold to make sparser solution
%     Iizscore = zscore(Ii(:,k));
%     Gizscore = zscore(Gi(:,k));
%     thresh = 0.5;
%     for i=1:length(I)
%         if (abs(Iizscore(i)) < thresh) && (abs(Gizscore(i)) < thresh)
%             Ii_clip(i,k) = (abs(Iizscore(i))/thresh)*Ii(i,k);
%             Gi_clip(i,k) = (abs(Gizscore(i))/thresh)*Gi(i,k);
%         else
%             Ii_clip(i,k) = Ii(i,k);
%             Gi_clip(i,k) = Gi(i,k);
%         end
%     end

    Iilfzscore = zscore(Iilf(:,k));
    Gilfzscore = zscore(Gilf(:,k));
    thresh = 0;
    for i=1:length(I)
        if (abs(Iilfzscore(i)) < thresh) && (abs(Gilfzscore(i)) < thresh)
            Iilf_clip(i,k) = (abs(Iilfzscore(i))/thresh)*Iilf(i,k);
            Gilf_clip(i,k) = (abs(Gilfzscore(i))/thresh)*Gilf(i,k);
        else
            Iilf_clip(i,k) = Iilf(i,k);
            Gilf_clip(i,k) = Gilf(i,k);
        end
    end
    clearvars W;
end


