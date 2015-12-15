clearvars -except traces segy cube Lcube_true

%matlabpool

tic;

i=5;
% from data fault cube
if i==0
    n_iline=segy{1}.n_iline;
    n_xline=segy{1}.n_xline;
    n_samples=500;  %segy{1}.n_samples;
    data=traces.data;
    
    
% model 1 - single noise-free elliptical plane
elseif i==1
    n_iline=100;
    n_xline=100;
    n_samples=100;
    % ellipse size and coordinates
    xc=50;
    yc=50;
    zc=50;
    a=50;
    b=50;
    % ellipse rotation
    theta=deg2rad(5);
    phi=deg2rad(0);
    PHI=deg2rad(90);
    [data,cube,Lcube_true]=model1(n_iline,n_xline,n_samples,xc,yc,zc,a,b,theta,phi,PHI);
    
    
% model 2 - single noisy elliptical plane
elseif i==2
    n_iline=100;
    n_xline=100;
    n_samples=100;    
    % ellipse size and coordinates
    xc=50;
    yc=50;
    zc=50;
    a=45;
    b=45;
    % ellipse rotation
    theta=deg2rad(120);
    phi=deg2rad(0);
    PHI=deg2rad(20);
    %fraction of noise introduced
    noise_frac=0.1;
    [data,cube,Lcube_true]=model2(n_iline,n_xline,n_samples,xc,yc,zc,a,b,theta,phi,PHI,noise_frac);   
    
    
% model 3 - multiple random, noise-free elliptical planes
elseif i==3
    n_iline=100;
    n_xline=100;
    n_samples=100;  
    %number of random faults generated
    nfaults=20;
    amax=10;
    bmax=10;
    [data,cube,Lcube_true]=model3(n_iline,n_xline,n_samples,nfaults,amax,bmax);

    
% model 4 - multiple random, noisy elliptical planes
elseif i==4
    n_iline=100;
    n_xline=100;
    n_samples=100;  
    nfaults=100;
    amax=50;
    bmax=50;
    noise_frac=0.1;
    [data,cube,Lcube_true]=model4(n_iline,n_xline,n_samples,nfaults,amax,bmax,noise_frac);
    

% data from cube
elseif i==5
    
    n_iline=size(cube,2);
    n_xline=size(cube,3);
    n_samples=size(cube,1);
    % unravel length cube into trace procession
    for i=1:n_iline
        ind=(1+(i-1)*n_xline):1:(i*n_xline);
        data(:,ind)=cube(:,i,:);
    end
end

% threshold segment length to consider in tracking
thres_len=3;

min_iline=1;
max_iline=n_iline;
min_xline=1;
max_xline=n_xline;

iLcube=zeros(n_samples,max_iline-min_iline+1,max_xline-min_xline+1);
xLcube=iLcube;

% get fault length sections for all inlines
for i=1:(max_iline-min_iline+1)
    ind=(min_xline+(min_iline+i-2)*n_xline):1:(max_xline+(min_iline+i-2)*n_xline);
    thin_slice=(data(1:n_samples,ind));
    iLcube(:,i,:)=IDcube_simple(thin_slice,thres_len);
end

% get fault length sections for all inlines
'inlines'
for i=min_iline:max_iline
    i
    ind=(min_xline+(i-1)*n_xline):1:(max_xline+(i-1)*n_xline);
    thin_slice=(data(1:n_samples,ind));
    iLcube(:,i-min_iline+1,:)=IDcube_simple(thin_slice,thres_len);
end

% get fault length sections for all crosslines
% parfor j=1:(max_xline-min_xline)
%     ind=(n_xline*(min_iline-1)+j+min_xline-1):n_xline:(n_xline*(max_iline-1)+j+min_xline-1);
%     thin_slice=(data(1:n_samples,ind));
%     xLcube(:,:,j)=IDcube_simple(thin_slice,thres_len);
% end

% get fault length sections for all crosslines
'crosslines'
for j=min_xline:max_xline
    j
    ind=(n_xline*(min_iline-1)+j):n_xline:(n_xline*(max_iline-1)+j);
    thin_slice=(data(1:n_samples,ind));
    xLcube(:,:,j-min_xline+1)=IDcube_simple(thin_slice,thres_len);
end

% match inline and crossline cubes and get longest lengths
ind1=iLcube>=xLcube;
ind2=xLcube>iLcube;
Lcube=ind1.*iLcube+ind2.*xLcube;

% go through time slices and group together faults using trace_wfit func
'slices'
Lcube2=zeros(size(Lcube));
for k=1:size(Lcube,1)
    k
    Lcube2(k,:,:)=IDslice(squeeze(Lcube(k,:,:)));
end

% unravel length cube into trace procession
% for i=min_iline:max_iline
%     ind=(min_xline+(i-1)*n_xline):1:(max_xline+(i-1)*n_xline);
%     Ltraces(:,ind)=Lcube(:,i-min_iline+1,:);
% end

toc

%matlabpool close