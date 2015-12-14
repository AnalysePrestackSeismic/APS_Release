function [] = quiverplot(data,v1,v2,subsample,fignum)
% plot data with quiver overlay
    [nZ,nX] = size(data);
    [xmesh,zmesh]=meshgrid((1:subsample:nX)',(1:subsample:nZ)');
    figure(fignum);
    imagesc(data); axis tight; colormap(gray);
    hold all;
    quiver(xmesh,zmesh,downsample(downsample(v1,subsample)',subsample)',downsample(downsample(v2,subsample)',subsample)',0.5,'r','MaxHeadSize',0.25)
    hold off
end