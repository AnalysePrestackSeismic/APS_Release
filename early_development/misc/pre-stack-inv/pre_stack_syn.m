% PreStack Inversion Synthetic

clear all

%loading file
W=load('BERG2b');

%assigning vectors
Depth=W(:,1);
Thick=W(:,2);
GR=W(:,3);
Sw=W(:,4);
Vp=W(:,5);
Vs=W(:,6);
RHOB=W(:,7);

angles = [5 10 15 20 25]; % angle stacks in degrees
angles = angles.*pi()/180;

n = length(Depth);

% Set up impedance matrix
imp(:,1) = Vp.*RHOB;
imp(:,2) = Vs.*RHOB;

del_imp = zeros(n,3);
for i=2:1:n
    del_imp(i,1) = (imp(i,1)-imp(i-1,1))/(imp(i,1)+imp(i-1,1)); % del_ip/2ip
    del_imp(i,2) = (imp(i,2)-imp(i-1,2))/(imp(i,2)+imp(i-1,2)); % del_is/2is
    del_imp(i,3) = (imp(i,2)+imp(i-1,2))/(imp(i,1)+imp(i-1,1)); % is/ip    
end

for i=1:1:length(angles)
    [wavelet(:,i)] = ricker_wavelet_function(40/i);
end

r = zeros(n,3);
for i=1:1:length(angles)
    r(:,i) = (1+tan(angles(i)).^2).*del_imp(:,1)-8.*(del_imp(:,3)).^2.*sin(angles(i)).^2.*del_imp(:,2);
    syn_obs(:,i) = conv2(r(:,i),wavelet(:,i),'same');
end

% Background model

prior{1} = fit(Depth,imp(:,1),'poly1');
prior{2} = fit(Depth,imp(:,2),'poly1');

prior_imp(:,1) = prior{1}.p1*Depth+prior{1}.p2;
prior_imp(:,2) = prior{2}.p1*Depth+prior{2}.p2;
prior_imp(:,3) = prior_imp(:,2)./prior_imp(:,1);

%
weight = [0.85 0.1 0.05];

j = vertcat(imp(:,1).*rand(1),imp(:,2).*rand(1));

loss = @(p)pre_stack_obj(syn_obs,prior_imp,weight,wavelet,angles,p');

[minimum,fval] = anneal(loss, j');

figure(1)
subplot(4,1,1); plot(imp(:,1)); hold on; plot(imp(:,2)); plot(prior_imp(:,1)); plot(prior_imp(:,2)); hold off;
subplot(4,1,2); imagesc(imp);
subplot(4,1,3); imagesc(syn_obs);
subplot(4,1,4); plot(minimum(1:length(imp(:,1)))); hold on; plot(minimum(length(imp(:,1))+1:end)); hold off;
