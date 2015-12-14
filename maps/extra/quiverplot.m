function [] = quiverplot(data,nX,nZ,v1,v2,l1,l2,subsample,fignum)
% plot data with quiver overlay
    [xmesh,zmesh]=meshgrid((1:subsample:nX)',(1:subsample:nZ)');
    figure(fignum);
    subplot(1,3,1); imagesc(data,[-1 1]); axis tight; colormap(gray);
    subplot(1,3,2); imagesc(data,[-1 1]); axis tight; colormap(gray);
    subplot(1,3,2); imagesc(l1-l2./(l1+l2)); axis tight;
    hold all;
    quiver(xmesh,zmesh,downsample(downsample(v1,subsample)',subsample)',downsample(downsample(v2,subsample)',subsample)',0.5,'r')
end