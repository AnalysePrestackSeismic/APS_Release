function [Iilf, Gilf, flag, relres, iter] = iginv_segy(seismic,Iprewhite,Gprewhite,Ismooth,Gsmooth)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nfiles = length(seismic);
npos = ((2*seismic{1}.aperture)+1)^2;
ntraces = seismic{1}.n_traces;
nt = seismic{1}.n_samples;
ntw = 101;
ntwh = floor((ntw+1)/2);
finish = 301;
Iilf = zeros(nt,finish);
Gilf = zeros(nt,finish);
angles = zeros(1,nfiles);

for ii = 1:nfiles
    angles(1,ii) = (sin(pi/180*seismic{ii}.angle)).*(sin(pi/180*seismic{ii}.angle));
end

read_idx(:,1:2) = [(1:npos:(ntraces*npos)-npos+1)',(npos:npos:ntraces*npos)'];

wrms = 1; % Scales the RMS of the wavelet to be equal to wrms
tic

for ij = 1:finish

    fprintf('Inverting trace %d of %d (%d%% complete).\n',ij,finish,round((ij/finish)*100));
    for jj = 1:nfiles
        tmp{jj} = segy_read_traces(seismic{jj}.filepaths,seismic{jj}.n_samples,seismic{jj}.process(read_idx(ij,1):read_idx(ij,2),:));
        S(:,read_idx(jj,1):read_idx(jj,2)) = tmp{jj}.data;
    end
    
    % Generate wavelets
    Sft_tmp = abs(fft(S))';
    Sft_tmp = reshape(Sft_tmp,npos,[]);
    Sft_tmp = mean(Sft_tmp,1);
    Sft_idx = logical(spdiags(ones(nfiles,nfiles*nt),(1:nfiles:nfiles*nt)-1,nfiles,nfiles*nt));
     
    for ik = 1:nfiles
        Sft(:,ik) = Sft_tmp(Sft_idx(ik,:))';
    end
      
    wraw = ifft(Sft/sqrt(nt),'symmetric');
        
    for kk = 1:nfiles
        wtmp = (2/sqrt(pi))*wraw(1:ntwh,kk).*cos(.5*pi*(0:ntwh-1)'/(ntwh-1));
        sn0 = wtmp(1)+2*sum(wtmp(2:end).*(-1).^(1:ntwh-1)');
        wtmp(1) = wtmp(1)-sn0;
        wtmp = [wtmp(ntwh:-1:2); wtmp];
        wtmp = wtmp/(norm(wtmp/wrms)/sqrt(length(wtmp)));
        if isnan(wtmp)
            w(:,kk) = 0;
        else
            w(:,kk) = wtmp;
        end
    end
    
    wangle = w*(spdiags(angles',0,nfiles,nfiles));
    
    % Whiten wavelets
    lambdaI = Iprewhite;
    lambdaG = Gprewhite;
    ww = w;
    ww(ntwh,:) = ww(ntwh,:)+lambdaI;
    wanglew = wangle;
    wanglew(ntwh,:) = wanglew(ntwh,:)+lambdaG;

    % Add linear trend to wavelets
    c = tan(18*pi/180);
    wlinfitw = ww*(spdiags((1-(angles'/c)),0,nfiles,nfiles));
    %wlinfitw = ww;
    
    for ji = 1:nfiles
        cols = [wlinfitw(:,ji); zeros(nt-1,1)];
        rows = zeros(nt,1);
        m = length(cols);
        x = [rows(nt:-1:2) ; cols(:)];                 
        cols_idx = (0:m-1)';
        rows_idx = nt:-1:1;
        
        Wlinfitw_tmp = cols_idx(:,ones(nt,1)) + rows_idx(ones(m,1),:);
        Wlinfitw_tmp(:) = x(Wlinfitw_tmp);
        Wlinfitw_tmp = Wlinfitw_tmp(((m-nt)/2)+1:m-((m-nt)/2),:);
        Wlinfitw(((npos*(ji-1))*nt)+1:npos*ji*nt,:) = repmat(Wlinfitw_tmp,npos,1);
        
        cols = [wanglew(:,ji); zeros(nt-1,1)];
        x = [rows(nt:-1:2) ; cols(:)];  
        Wanglew_tmp = cols_idx(:,ones(nt,1)) + rows_idx(ones(m,1),:);
        Wanglew_tmp(:) = x(Wanglew_tmp);
        Wanglew_tmp = Wanglew_tmp(((m-nt)/2)+1:m-((m-nt)/2),:);
        Wanglew(((npos*(ji-1))*nt)+1:npos*ji*nt,:) = repmat(Wanglew_tmp,npos,1);
    end
  
    % Tikhonov Regularization Matrix
    B = [-0.5,1,-0.5];
    B = repmat(B,2*nt,1);
    L = chol(spdiags(B,[-1,0,1],2*nt,2*nt));
    Iscale = Ismooth;
    Gscale = Gsmooth;
    L(1:nt,:) = Iscale*L(1:nt,:);
    L(nt+1:end,:) = Gscale*L(nt+1:end,:);
    op_linfit = [Wlinfitw,Wanglew;L];
    
    [inv_linfit,flag(ij),relres(ij),iter(ij)] = lsqr(op_linfit,[reshape(S,[],1);zeros(2*nt,1)],1e-2,100,[],[],[S(:,(npos+1)/2);zeros(nt,1)]);
    %[inv_linfit] = lsqr(op_linfit,[reshape(S,[],1);zeros(2*nt,1)],1e-2,100,[],[],[S(:,(npos+1)/2);-S(:,(npos+1)/2)/tan(20*pi/180)]);
    
    Iilf(:,ij) = inv_linfit(1:nt);
    Gilf(:,ij) = inv_linfit(nt+1:end)-(Iilf(:,ij)/c);
        
    fprintf('Elapsed time: %.1f minutes.\n',toc/60);
    fprintf('Time per trace: %.2f seconds.\n',toc/ij);
    fprintf('Time to completion: %.1f minutes.\n\n',((toc/ij)*(finish-ij))/60);
end

end

