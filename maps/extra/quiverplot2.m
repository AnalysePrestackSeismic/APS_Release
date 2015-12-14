function [] = quiverplot2(data,nX,nZ,v1,v2,v3,v4,l1,l2,subsample,fignum)
% plot data with quiver overlay
    [xmesh,zmesh]=meshgrid((1:subsample:nX)',(1:subsample:nZ)');
    figure(fignum);
    subplot(3,2,1); imagesc(data,[-0.1 0.1]); axis tight; colormap(gray);

    subplot(3,2,2); imagesc(l1-l2./(l1+l2)); axis tight; colormap(gray);

    subplot(3,2,3); imagesc(data,[-0.1 0.1]); axis tight; colormap(gray);
    
    hold all;
    quiver(xmesh,zmesh,downsample(downsample(v1,subsample)',subsample)',downsample(downsample(v2,subsample)',subsample)',0.5,'r')
    hold off; 
    
    subplot(3,2,4); imagesc(atand(v2./v1),[-45 45]); axis tight; colormap(gray);    
    
    subplot(3,2,5); imagesc(data,[-0.1 0.1]); axis tight; colormap(gray);
    
    hold all;
    quiver(xmesh,zmesh,downsample(downsample(v3,subsample)',subsample)',downsample(downsample(v4,subsample)',subsample)',0.5,'r')
    hold off; 
end