function [x] = annealing(angle_stacks,prior,weight,wavelet,angles,p)

w = [1 0 0]; % weighting in objective function terms

% Upper and Lower bounds on Ip and Is
upper = ;

lower = ;



%pre_stack_obj(syn_obs,prior,w,wavelet,angles,imp)
loss = @(p)pre_stack_obj(syn_obs,prior,w,wavelet,angles,p);
guess(:,1) = 2000.*ones(n,1);
guess(:,2) = 1000.*ones(n,1);
guess = [guess(:,1) guess(:,2)];
[x fval] = simulannealbnd(loss,guess)
imp

end