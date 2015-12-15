function f = pre_stack_obj(angle_stacks,prior,weight,wavelet,angles,imp)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% 
% Inputs:
%   angle_stacks   weights
%   synthetic
% Impedance matrix (imp) - Col1 = Ip, Col2 = Is

[n_samples,n_angles] = size(angle_stacks);
[~,n_wavelets] = size(wavelet);
% [l_wavelets,n_wavelets] = size(wavelets);


% term 1
% generate synthetics
del_imp = zeros(n_samples,n_angles);

for i=2:1:n_samples
    del_imp(i,1) = (imp(i,1)-imp(i-1,1))/(imp(i,1)+imp(i-1,1)); % del_ip/2ip
    del_imp(i,2) = (imp(i,2)-imp(i-1,2))/(imp(i,2)+imp(i-1,2)); % del_is/2is
    del_imp(i,3) = (imp(i,2)+imp(i-1,2))/(imp(i,1)+imp(i-1,1)); % is/ip    
end

r = zeros(n_samples,n_angles);
model_stacks = zeros(n_samples,n_angles);
for i=1:1:n_angles
    r(:,i) = (1+tan(angles(i)).^2).*del_imp(:,1)-8.*(del_imp(:,3)).^2.*sin(angles(i)).^2.*del_imp(:,2);
    model_stacks(:,i) = conv(r(:,i),wavelet(:,i),'same');
end

term1_top = abs(angle_stacks-model_stacks);
term1_top = sum(sum(term1_top,1),2);

term1_bot = abs(angle_stacks);
term1_bot = sum(sum(term1_bot,1),2);

term1 = (term1_top/term1_bot);

% term2
% prior: Col1 = Ip, Col2 = Is
term2a_top = abs(prior(:,1)-imp(:,1));
term2a_top = sum(term2a_top,1);
term2a_bot = abs(prior(:,1));
term2a_bot = sum(term2a_bot,1);

term2b_top = abs(prior(:,2)-imp(:,2));
term2b_top = sum(term2b_top,1);
term2b_bot = abs(prior(:,2));
term2b_bot = sum(term2b_bot,1);

term2 = (term2a_top/term2a_bot)+(term2b_top/term2b_bot);

% term3
% prior: Col3 = Is/Ip
term3_top = abs(prior(:,3)-del_imp(:,3));
term3_top = sum(term3_top,1);

term3_bot = abs(prior(:,3));
term3_bot = sum(term3_bot,1);

term3 = term3_top/term3_bot;

% objective function
f = weight(1)*term1+weight(2)*term2+weight(3)*term3;

end

