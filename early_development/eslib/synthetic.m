function [syn] = synthetic(ref,wavelet)

[t,x]=size(ref);
for i=1:x
    syn=conv2(ref(:,i),wavelet,'same');
end
