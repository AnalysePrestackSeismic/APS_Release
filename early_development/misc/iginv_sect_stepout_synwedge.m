% Some setup parameters
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

nt = 64;
taper = 5; % Defines cos-taper length in samples for beginning and end of trace (might not be needed)
wrms = 1; % Scales the RMS of the wavelet to be equal to wrms
stepout = 3; % Trace stepout for inversion (total traces used equals 1+2*stepout)
ntrace = 1+2*stepout;
residual = nan(length(near)-2*stepout);
iter = nan(length(near)-2*stepout);
chi = nan(length(near)-2*stepout);

% Estimate wavelets for near, mid and far
% nearft = mean(abs(fft(near)),2);
% wraw = ifft(nearft/sqrt(nt),'symmetric');
% wtmp = (2/sqrt(pi))*wraw(1:51).*cos(.5*pi*(0:50)'/(50));
% sn0 = wtmp(1)+2*sum(wtmp(2:end).*(-1).^(1:50)');
% wtmp(1) = wtmp(1)-sn0;
% wtmp = [wtmp(51:-1:2); wtmp];
% wtmp = wtmp/(norm(wtmp/wrms)/sqrt(length(wtmp)));
% wnear = wtmp;
% 
% midft = mean(abs(fft(mid)),2);
% wraw = ifft(midft/sqrt(nt),'symmetric');
% wtmp = (2/sqrt(pi))*wraw(1:51).*cos(.5*pi*(0:50)'/(50));
% sn0 = wtmp(1)+2*sum(wtmp(2:end).*(-1).^(1:50)');
% wtmp(1) = wtmp(1)-sn0;
% wtmp = [wtmp(51:-1:2); wtmp];
% wtmp = wtmp/(norm(wtmp/wrms)/sqrt(length(wtmp)));
% wmid = wtmp;
% 
% farft = mean(abs(fft(far)),2);
% wraw = ifft(farft/sqrt(nt),'symmetric');
% wtmp = (2/sqrt(pi))*wraw(1:51).*cos(.5*pi*(0:50)'/(50));
% sn0 = wtmp(1)+2*sum(wtmp(2:end).*(-1).^(1:50)');
% wtmp(1) = wtmp(1)-sn0;
% wtmp = [wtmp(51:-1:2); wtmp];
% wtmp = wtmp/(norm(wtmp/wrms)/sqrt(length(wtmp)));
% wfar = wtmp;

% w = [wnear,wmid,wfar];

% Make angles
angles = [5,15,25,35];
nangle = length(angles);
angles_stepout = sort(repmat(angles,1,ntrace));
sin2theta = (sin(angles.*pi/180)).*(sin(angles.*pi/180));
sin2theta_stepout = (sin(angles_stepout.*pi/180)).*(sin(angles_stepout.*pi/180));

% Start of inversion loop
for k = 1+stepout:length(near)-stepout
%for k = 12:12
    fprintf('\nInverting traces %d to %d of %d\n',k-stepout,k+stepout,length(near))
    
    % Taper start and end of traces
    near_taper = near(:,k-stepout:k+stepout);
    mid_taper = mid(:,k-stepout:k+stepout);
    far_taper = far(:,k-stepout:k+stepout);
    ufar_taper = ufar(:,k-stepout:k+stepout);
    
    for i = 1:(2*stepout)+1
        near_taper(end-taper:end,i) = (2/sqrt(pi))*near(end-taper:end,k-stepout+i-1).*cos(.5*pi*(0:taper)'/taper);
        near_taper(1:taper+1,i) = (2/sqrt(pi))*near(1:taper+1,k-stepout+i-1).*cos(.5*pi*(taper:-1:0)'/taper);
        mid_taper(end-taper:end,i) = (2/sqrt(pi))*mid(end-taper:end,k-stepout+i-1).*cos(.5*pi*(0:taper)'/taper);
        mid_taper(1:taper+1,i) = (2/sqrt(pi))*mid(1:taper+1,k-stepout+i-1).*cos(.5*pi*(taper:-1:0)'/taper);
        far_taper(end-taper:end,i) = (2/sqrt(pi))*far(end-taper:end,k-stepout+i-1).*cos(.5*pi*(0:taper)'/taper);
        far_taper(1:taper+1,i) = (2/sqrt(pi))*far(1:taper+1,k-stepout+i-1).*cos(.5*pi*(taper:-1:0)'/taper);
        ufar_taper(end-taper:end,i) = (2/sqrt(pi))*ufar(end-taper:end,k-stepout+i-1).*cos(.5*pi*(0:taper)'/taper);
        ufar_taper(1:taper+1,i) = (2/sqrt(pi))*ufar(1:taper+1,k-stepout+i-1).*cos(.5*pi*(taper:-1:0)'/taper);
    end
    
    S = [near_taper,mid_taper,far_taper,ufar_taper];

    % Estimate wavelets
    Sft = abs(fft(S));
    wraw = ifft(Sft/sqrt(nt),'symmetric');
    for i=1:ntrace*nangle
        wtmp = (2/sqrt(pi))*wraw(1:51,i).*cos(.5*pi*(0:50)'/(50.5));
        sn0 = wtmp(1)+2*sum(wtmp(2:end).*(-1).^(1:50)');
        wtmp(1) = wtmp(1)-sn0;
        wtmp = [wtmp(51:-1:2); wtmp];
        wtmp = wtmp/(norm(wtmp/wrms)/sqrt(length(wtmp)));
        w(:,i) = wtmp;
    end
    
%     Sft = abs(fft(S));
%     wraw = ifft(Sft/sqrt(nt),'symmetric');
%     for i=1:nangle
%         wtmp = (2/sqrt(pi))*wraw(1:51,i).*cos(.5*pi*(0:50)'/(50.5));
%         sn0 = wtmp(1)+2*sum(wtmp(2:end).*(-1).^(1:50)');
%         wtmp(1) = wtmp(1)-sn0;
%         wtmp = [wtmp(51:-1:2); wtmp];
%         wtmp = wtmp/(norm(wtmp/wrms)/sqrt(length(wtmp)));
%         w(:,i) = wtmp;
%     end

    % Make reflection coeffs
    % for i=1:nangle
    %     R(:,i) = I+G*sin2theta(i);
    % end

    % Make true synthetic
    % for i=1:nangle
    %     S(:,i) = conv2(R(:,i),w(:,i),'same');
    % end

    % Make 'conditioned' synthetic
    % for i=1:nangle
    %     Sc(:,i) = conv2(R(:,i),w(:,2),'same');
    % end

    % Estiamte I and G from line fitting to the data
%     for i=1:nt
%         tmp = robustfit(sin2theta_stepout,S(i,:),'ols');
%         Is(i,k) = tmp(1,1);
%         Gs(i,k) = tmp(2,1);
%     end
    
    for i=1:nt
        tmp = robustfit(sin2theta,[S(i,stepout+1),S(i,1+2*(stepout+1)),S(i,1+3*(stepout+1)),S(i,1+2*(stepout+1))],'ols');
        Is(i,k) = tmp(1,1);
        Gs(i,k) = tmp(2,1);
    end
    
%     % Estiamte I and G from line fitting to the data
%     for i=(stepout*nt)+1:(stepout+1)*nt
%         tmp = robustfit(sin2theta,S(i,:),'ols');
%         Is(i-stepout*nt,k) = tmp(1,1);
%         Gs(i-stepout*nt,k) = tmp(2,1);
%     end

    % Estimate wavelet scalers
    % for i=1:nangle
        %scale(1,i) = (norm(S)/(norm(R)*norm(w)));
        %scale = [1,1,1];
    % end

    % Make wavelet convolution matrix
    for i = 1:ntrace*nangle
        W(:,nt*(i-1)+1:nt*i) = convmtx(w(:,i),nt);
    end
    clip = (size(W,1)-nt)/2;
    W = W(clip+1:end-clip,:);
    W = sparse(W);

    % Make whitenning matrix
    tmp_white = eye(nt);
    white = repmat(tmp_white,ntrace*nangle,2);
    lambdaI = 1;
    lambdaG = 1;
    white(:,1:nt) = white(:,1:nt)*lambdaI;
    white(:,1+nt:end) = white(:,1+nt:end)*lambdaG;
    white = sparse(white);

    % Make operator without linear trend
    op = [W',W'];
    for i=1:ntrace*nangle
        op(1+(i-1)*nt:i*nt,1+nt:end) = op(1+(i-1)*nt:i*nt,1+nt:end)*sin2theta_stepout(i);
    end
    op = op+white;

    % Estimate c
    minchi = 10;
    maxchi = 60;
    carray = (tan((minchi:5:maxchi)*pi/180))';
    relres = zeros(1+(maxchi/5-minchi/5),1);
    for n=1:(maxchi/5-minchi/5)+1
        op_linfit = op;
        for i=1:ntrace*nangle
            op_linfit(1+(i-1)*nt:i*nt,1:nt) = op_linfit(1+(i-1)*nt:i*nt,1:nt)*(1-(sin2theta_stepout(i)/carray(n)));
        end
        %op_linfit = op_linfit+white;

        % Do inversion to find best c
        [~,~,relres(n,1)] = lsqr(op_linfit,reshape(S,[],1),1e-3,5,[],[],[Is(:,k);zeros(nt,1)]);
        %relrestrack(n,1) = relres;
        %relRMSdG(n,1) = (norm(inv_linfit(1+length(S):end))/sqrt(length(S)))/(norm(S)/sqrt(length(S)*nangle));
    end
    relresi = interpft(relres,51);
    [~,chi(k,1)] = min(relresi);
    fprintf('chi = %d\n',chi(k)+minchi-1);
    c = tan((chi(k)+minchi-1)*pi/180);
%     c = tan(16*pi/180);
    
    % Make operater with linear trend
    op_linfit = op;
    for i=1:ntrace*nangle
        op_linfit(1+(i-1)*nt:i*nt,1:nt) = op_linfit(1+(i-1)*nt:i*nt,1:nt)*(1-(sin2theta_stepout(i)/c));
    end
    
    % Tikhonov Regularization Matrix
    B = [-0.5,1,-0.5];
    B=repmat(B,2*251,1);
    %scale = (1:2*nt)'/5;
    %B = B.*repmat(scale,1,3);
    %[~,L] = lu(spdiags(B,[-1,0,1],2*nt,2*nt));
    L = chol(spdiags(B,[-1,0,1],2*nt,2*nt));
    Iscale = 40;
    Gscale = 10;
    L(1:nt,:) = Iscale*L(1:nt,:);
    L(nt+1:end,:) = Gscale*L(nt+1:end,:);
    op_linfit = [op_linfit;L];
    
    % Do inversions
    % inv = lsqr(op,reshape(S,[],1),1e-6,500,[],[],[S(:,1);-S(:,1)/c]);
    [inv_linfit,~,residual(k,1),iter(k,1)] = lsqr(op_linfit,[reshape(S,[],1);Is(:,k);zeros(nt,1)],1e-3,1000,[],[],[Is(:,k);zeros(nt,1)]);
    fprintf('convergence after %d iterations with relative residual of %f\n',iter(k),residual(k));

%     Ii(:,k) = inv(1:nt,1);
%     Gi(:,k) = inv(1+nt:end,1);

    % Select the kth trace as the result for the kth position
    Iilf(:,k) = inv_linfit(1:nt);
    Gilf(:,k) = inv_linfit(nt+1:end)-(Iilf(:,k)/c);

    % Average the kth +/- stepout traces for the results at the kth position
%     Iilf(:,k) = mean(reshape(inv_linfit(1:length(S),1),[],ntrace),2);
%     Gilf(:,k) = mean(reshape(inv_linfit(1+length(S):end,1),[],ntrace),2)-(Iilf(:,k)/c);

    % Threshold to make sparser solution
%     Iizscore = zscore(Ii(:,k));
%     Gizscore = zscore(Gi(:,k));
%     thresh = 0.5;
%     for i=1:nt
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
    for i=1:nt
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


