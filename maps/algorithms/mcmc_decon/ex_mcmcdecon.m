n = 1024;
%dt = 0.002;
%fsample = 1./dt;
rand('state',0)
randn('state',3)

% create true wavelet

wr = bliss_ricker(n,1/1000,60,45)' + bliss_ricker(n,1/1000,120,0)';
w = [wr(end-19:end) wr(1:21)];
m = 41;					% wavelet length

% create replectivity

np = n + m-1;
p = input(['Choice of reflectivity distribution\n',...
	   '0 < p < 1 : Bernoulli-Gaussian with probability p\n',...
	   'p = 1 : bilateral exponential\n',...
	   'p = 2, 3, ...: sign(z)abs(z)^p with z normal\n',...
	   'your choice ? ']);
if p > 1
  r = randn(1,np);
  r = sign(r).*abs(r).^p;
elseif p == 1
  r = 2*rand(1,np) - 1;
  r = sign(r).*log(1 - abs(r));
elseif p > 0
  tmp = find(rand(1,np) < p);
  r = zeros(1,np);
  r(tmp) = randn(1,length(tmp));
else
  fprintf('wrong choice\n')
  return
end
r = r./sqrt((r*r')/np);			% nomalise to have power 1

% set noise level
snr = input('S/N (standard dev. of signal/standard dev. noise) ratio ? ');
sn = sqrt((w*w'))/snr			% noise staandard deviation;

% create data
data = filter(w,1,r);
data = data(m:np) + sn*randn(1,n);

%%% deconvolution %%%

nc = input('Numbers of warmup and simulation cycles (in bracket []) ');

[rhat, rm] = mcmcdecon(data, w, sn, nc(1), nc(2), true);

c0 = (r*rhat)./sqrt((rhat'*rhat)*(r*r'));
c1 = (r*rm)./sqrt((rm'*rm)*(r*r'));
fprintf('corr. between decon output  & reflect. %.4g\n', c1)
fprintf('corr. between Wiener output & reflect. %.4g\n', c0)

figure(1)
subplot(4,1,1); plot(r); axis tight; title('true reflectivity')
subplot(4,1,2); plot(data); axis tight; title('observation')
subplot(4,1,3); plot(rhat); axis tight;
title(sprintf('Wiener decon., corr. = %.4g', c0))
subplot(4,1,4); plot(rm); axis tight;
title(sprintf('MCMC decon., corr. = %.4g', c1))

figure(2); plot(1:np,r, 1:np,rm); axis tight
legend('true', 'decon.')
