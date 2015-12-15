function [tens] = struc_tens_2d(data)
% calculate structural tensors from a 2d image. 

    han0 = waitbar(0,'Initializing...','Name','Computing structural tensors');
    [nZ,nX] = size(data);
    
    % calculate directional derivatives using Gaussian derivative filters
    [Ix,Iz] = gradient(data);

    % make smoothed structure tensor components
    H = makeGaussianFilter(7);
    Ixz =  reshape(conv2(H',H,Ix.*Iz,'same'),[],1);
    Ixx =  reshape(conv2(H',H,Ix.*Ix,'same'),[],1);
    Izz =  reshape(conv2(H',H,Iz.*Iz,'same'),[],1);
    
    % do eigenvalue decomposition
    tic
    for ii = 1:length(Izz)
        [e1{ii},e2{ii},l1(ii),l2(ii)] = eigen_decomposition([Ixx(ii),Ixz(ii);Ixz(ii),Izz(ii)]);
        waitbar(ii/length(Izz),han0,sprintf('%.1f%% complete, approx. %.0f minutes remaining',100*ii/length(Izz),((toc/60)/ii)*(length(Izz)-ii)));
    end
    toc;
    close(han0);

    % eigenvalues (l1>l2)
    l1 = reshape(l1,nZ,nX);
    l2 = reshape(l2,nZ,nX);
    
    % x and z components of eigenvectors for l1
    tmp = cell2mat(e1);
    e1x = reshape(tmp(1,:),nZ,nX);
    e1z = reshape(tmp(2,:),nZ,nX);
    
    % x and z components of eigenvectors for l2
    tmp = cell2mat(e2);
    e2x = reshape(tmp(1,:),nZ,nX);
    e2z = reshape(tmp(2,:),nZ,nX);
    
    % 3d matrix of eigenvalues and eigenvectors
    tens = cat(3,l1,e1x,e1z,l2,e2x,e2z);
    
    % make some plots
    quiverplot(data,nX,nZ,e2x,e2z,5,1); % plot data with quiver overlay
    
    data_smoothed = reshape(pcg(@(x)smooth(x,1000,nX,nZ,sign(e1x).*e1z),data(:),[],500),nZ,nX); % smooth along structure
    figure(2); subplot(1,3,1); imagesc(data_smoothed,[-1 1]); axis equal; axis tight; colormap(gray);
    similarity = sim(data,nX,nZ,tens);
    data_smoothed_with_sim = reshape(pcg(@(x)smooth(x,1000*(similarity.^20),nX,nZ,sign(e1x).*e1z),data(:),[],500),nZ,nX); % smooth along structure, weighted by similarity raised to a big power (e.g. Hale uses 8 in one paper)
    figure(2); subplot(1,3,2); imagesc(data_smoothed_with_sim,[-1 1]); axis equal; axis tight; colormap(gray);
    figure(2); subplot(1,3,3); imagesc(data_smoothed-data_smoothed_with_sim,[-0.2 0.2]); axis equal; axis tight; colormap(gray); % differece between smoothing with and without similarity
    
    figure(3); imagesc(atand(e2z./e2x),[-45 45]); axis equal; axis tight; % plot dip from structural tensors
end

function [H] = makeGaussianFilter(sigma)
% calculate normalised Gaussian filter
    siz = sigma*6;
    x=(-ceil(siz/2):ceil(siz/2))';
    H = exp(-(x.^2/(2*sigma^2)));
    H = H/sum(H(:));
end

function [e1,e2,l1,l2] = eigen_decomposition(M)
% calculate eigenvalues and eigenvectors
    [V D] = eig(M);
    [l,idx] = sort(abs(diag(D,0)),'descend');
    l1 = l(1);
    l2 = l(2);
    e1 = V(:,idx(1));
    e2 = V(:,idx(2));
end

function [] = quiverplot(data,nX,nZ,v1,v2,subsample,fignum)
% plot data with quiver overlay
    [xmesh,zmesh]=meshgrid((1:subsample:nX)',(1:subsample:nZ)');
    figure(fignum);
    imagesc(data,[-1 1]); axis equal; axis tight; colormap(gray);
    hold all;
    quiver(xmesh,zmesh,downsample(downsample(v1,subsample)',subsample)',downsample(downsample(v2,subsample)',subsample)',0.5,'r','MaxHeadSize',0.25)
    hold off
end

function [g] = tensfiltA(data,nX,nZ,u2)
% apply A'A filter from D. Hale, 2007, CWP Report 567, "Local dip filtering
% with directional Laplacians"
    for ii = 2:nZ
        for jj=2:nX
            u2i = u2(ii,jj);
            u1i = sqrt(1-u2i*u2i);
            a11 = 1-u1i*u1i;
            a12 = -u1i*u2i;
            a22 = 1-u2i*u2i;
            fa = data(ii,jj)-data(ii-1,jj-1);
            fb = data(ii,jj-1)-data(ii-1,jj);
            f1 = 0.5*(fa-fb);
            f2 = 0.5*(fa+fb);
            g1 = a11*f1+a12*f2;
            g2 = a12*f1+a22*f2;
            ga = 0.5*(g1+g2);
            gb = 0.5*(g1-g2);
            g(ii,jj) = ga;
            g(ii-1,jj-1) = g(ii-1,jj-1)-ga;
            g(ii,jj-1) = g(ii,jj-1)-gb;
            g(ii-1,jj) = g(ii-1,jj)+gb;
        end
    end
end

function [y] = smooth(x,alpha,nX,nZ,u2)
% smoothing by inversion regularization using Claerbout's wavekill filter
% (filter A from D. Hale, 2007, CWP Report 567, "Local dip filtering
% with directional Laplacians"). Alpha controls smoothness.
    x = reshape(x,nZ,nX);
    y = reshape(x+alpha.*tensfiltA(x,nX,nZ,u2),[],1);
end

function [y] = sim(data,nX,nZ,tens) 
% structure-orientated similarity from D. Hale, 2009, CWP Report 635,
% "Structure-oriented smoothing and semblance").
    y1 = pcg(@(x)smooth(x,3,nX,nZ,sign(tens(:,:,2)).*tens(:,:,3)),data(:),[],500);
    y1 = reshape(pcg(@(x)smooth(x,45,nX,nZ,sign(tens(:,:,5)).*tens(:,:,6)),y1.^2,[],500),nZ,nX); 
    y2 = pcg(@(x)smooth(x,3,nX,nZ,sign(tens(:,:,2)).*tens(:,:,3)),data(:).^2,[],500);
    y2 = reshape(pcg(@(x)smooth(x,45,nX,nZ,sign(tens(:,:,5)).*tens(:,:,6)),y2,[],500),nZ,nX);
    y = y1./y2;
    y(y<0) = 0;
    y(y>1) = 1;
end