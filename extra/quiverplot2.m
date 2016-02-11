function [] = quiverplot2(data,nX,nZ,v1,v2,v3,v4,l1,l2,subsample,fignum)
%% ------------------ Disclaimer  ------------------
% 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) makes no representation or warranty, express or implied, in 
% respect to the quality, accuracy or usefulness of this repository. The code
% is this repository is supplied with the explicit understanding and 
% agreement of recipient that any action taken or expenditure made by 
% recipient based on its examination, evaluation, interpretation or use is 
% at its own risk and responsibility.
% 
% No representation or warranty, express or implied, is or will be made in 
% relation to the accuracy or completeness of the information in this 
% repository and no responsibility or liability is or will be accepted by 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) in relation to it.
%% ------------------ License  ------------------ 
% GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
%% github
% https://github.com/AnalysePrestackSeismic/
%% ------------------ FUNCTION DEFINITION ---------------------------------
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