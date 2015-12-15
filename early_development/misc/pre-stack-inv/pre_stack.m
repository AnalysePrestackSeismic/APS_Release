% Prestack inversion via simulated annealing (Ma, 2002)

%

clear all

% Generate synthetic data
n = 1; % number of samples
dt = 0.004; % sample rate in s
angles = [5 20 30]; % angle stacks in degrees
angles = angles.*pi()/180;

% generate wavelets
% Set effective support and grid parameters.
lb = -8; ub = 8; w_l = 1;

% Compute Gaussian wavelet of order 8.
[psi,x] = gauswavf(lb,ub,w_l,2);
for i=1:1:length(angles)
    wavelet(:,i) = psi'.*(1/i);
end

%figure(1)
%plot(psi)

% Impedance matrix - Col1 = Ip, Col2 = Is
ip_is = 2;
imp = zeros(n,2);
% random impedances
% rng('default')
ip_min = 5000;
ip_max = 9000;
is_min = ip_min/ip_is;
is_max = ip_max/ip_is;
for i=1:1:n
    imp(i,1) = (ip_min + (ip_max-ip_min).*rand(1));
    imp(i,2) = (is_min + (is_max-is_min).*rand(1));
end

% for i=1:1:n
%     imp(i,1) = ip_min;
%     imp(i,2) = is_min;
% end
% 
% for i=250:1:n
%     imp(i,1) = ip_max;
%     imp(i,2) = is_max;
% end

del_imp = zeros(n,3);
for i=2:1:n
    del_imp(i,1) = (imp(i,1)-imp(i-1,1))/(imp(i,1)+imp(i-1,1)); % del_ip/2ip
    del_imp(i,2) = (imp(i,2)-imp(i-1,2))/(imp(i,2)+imp(i-1,2)); % del_is/2is
    del_imp(i,3) = (imp(i,2)+imp(i-1,2))/(imp(i,1)+imp(i-1,1)); % is/ip    
end

r = zeros(n,3);
for i=1:1:length(angles)
    r(:,i) = (1+tan(angles(i)).^2).*del_imp(:,1)-8.*(del_imp(:,3)).^2.*sin(angles(i)).^2.*del_imp(:,2);
    syn_obs(:,i) = conv2(r(:,i),wavelet(:,i),'same');
end

w = [0.7 0.2 0.1]; % weighting in objective function terms

prior(:,1) = ones(n,1)*mean(imp(:,1));
prior(:,2) = ones(n,1)*mean(imp(:,2));
prior(:,3) = prior(:,2)./prior(:,1);

% Upper and Lower bounds on Ip and Is
% prior: Col1 = Ip, Col2 = Is
upper(:,1) = prior(:,1)+3000; % Col1: Ip, Col2: Is
upper(:,2) = prior(:,2)+2000; 

lower(:,1) = prior(:,1)-3000; % Col1: Ip, Col2: Is
lower(:,2) = prior(:,2)-2000;

% initial guess
% guess(:,1) = 5000*ones(n,1);
% guess(:,2) = 3000*ones(n,1);
% 
% cost = @(p)pre_stack_obj(syn_obs,prior,w,wavelet,angles,p);
% 
% current_state = guess;
% current_cost = cost(guess);
% temp = 1;
% 
% m = 1;
% while m<1000;
%   
%     % Initalise a new state
%     state = zeros(n,2);
%     for i=1:1:n
%         state(i,1) = lower(i,1) + (upper(i,1)-lower(i,1)).*rand(1);
%         state(i,2) = lower(i,2) + (upper(i,2)-lower(i,2)).*rand(1);
%     end
%         
%     new_cost = cost(state);
%     
%     if (current_cost-new_cost) <= 0
%         current_state = state;
%     elseif (current_cost-new_cost)/temp > rand(1)
%         current_state = state;    
%     end
%     
%     %m_cost(m) = cost(state);
%     
%     temp = 0.8*temp;
%     
%     m = m + 1;
%         
% end
% % 
% figure(1)
% subplot(2,2,1)
% plot(imp(:,1));
% hold on
% plot(upper(:,1));
% plot(lower(:,1));
% hold off
% 
% subplot(2,2,2)
% plot(imp(:,2));
% hold on
% plot(upper(:,2));
% plot(lower(:,2));
% hold off
% 
% subplot(2,2,3)
% plot(imp);
% 
% subplot(2,2,4)
% plot(state);

%pre_stack_obj(syn_obs,prior,w,wavelet,angles,imp)
loss = @(p)pre_stack_obj(syn_obs,prior,w,wavelet,angles,p);
k = 1;
while k < 10000
    for i=1:1:n
        state(i,1) = lower(i,1) + (upper(i,1)-lower(i,1)).*rand(1);
        state(i,2) = lower(i,2) + (upper(i,2)-lower(i,2)).*rand(1);
    end
    cost(k) = loss(state);
    g(k,:) = state;
    k = k + 1;
end
%guess(:,1) = 3000.*ones(n,1);
%guess(:,2) = 2000.*ones(n,1);
%guess = [guess(:,1) guess(:,2)];
%[x fval] = simulannealbnd(loss,guess)
%imp

%[x] = SIMPSA('loss',guess,lower,upper);

%end

%figure(2)
%imagesc(syn_obs)

% loss = @(p)camel(p(1),p(2));
%loss = @(p)pre_stack_obj(syn_obs,prior,w,wavelet,angles,p);

%k = [imp(:,1) ; imp(:,2)]
%for i=1:1:500;
%    res(i) = loss(250*k/i);
%end

%[x f] = anneal(loss,[7000 8000 2850 1000]);
%value = loss(x)
%x'
%k = [imp(:,1) ; imp(:,2)]

% objective functio
% Import angle stacks
% Import wavelets

