% Testing the nonlinear fit tool: nlinfit
%

clear all

% Our model: Vexp(t) = Va*Vinf/(Va+(Vinf-Va)*exp(-ka*t*Vinf/(Vinf-Va)))
V_inf = 6000;
mdl = @(a,t)((a(1)*V_inf)./(a(1)+(V_inf-a(1))*exp((-a(2)*t*V_inf)/(V_inf-a(1)))));

% make some synthetic data
rng(9845,'twister')

% make a few different models
m(:,1) = [2000 ; 0.5];
m(:,2) = [2010 ; 0.52];
m(:,3) = [2050 ; 0.54];
% m(:,4) = [1990 ; 0.49];

t = [0:0.25:5]; %exprnd(2,100,1);
epsn = normrnd(0,10,length(t),1);

V_exp(:,1) = mdl(m(:,1),t) + epsn';
V_exp(:,2) = mdl(m(:,2),t) + epsn';
V_exp(:,3) = mdl(m(:,3),t) + epsn';
% V_exp(:,4) = mdl(m(:,4),t) + epsn';
%
V_exp(:,1) = V_exp(:,1);
V_exp(:,2) = V_exp(:,2);
V_exp(:,3) = V_exp(:,3);
% V_exp(:,4) = 

a0 = [1200; 0.3];

V_exp_inv = V_exp(1:end);

[ahat(:,1),r,J,cov,mse] = nlinfit(t,V_exp(:,1)',mdl,a0);
[ahat(:,2),r,J,cov,mse] = nlinfit(t,V_exp(:,2)',mdl,a0);
% ahat
% % 
trange = min(t):.01:max(t);
hold on
scatter(t,V_exp(:,1))
scatter(t,V_exp(:,2))
plot(trange,mdl(ahat(:,1),trange),'r')
plot(trange,mdl(ahat(:,2),trange),'r')
hold off
 
% % using weights
% 
%w = [1 1 1 1 1 1 1 1 1 1 1];
%w = w/mean(w);
%V_exp_w = sqrt(w).*V_exp;
% modelFunw = @(a,t) sqrt(w).*mdl(a,t);
% [bFitw,rw,Jw,Sigmaw,msew] = nlinfit(t,V_exp_w,modelFunw,a0);
% bFitw
% 
% figure
% trange = min(t):.01:max(t);
% hold on
% scatter(t,V_exp_w)
% plot(trange,mdl(bFitw,trange),'r')
% hold off