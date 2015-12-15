function synthetic = gen_smod(IpMod,IsMod,wavelet,angle)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%del_imp = NaN(length(IpMod),3);

for i=2:1:length(IpMod)
    del_ip(i,1) = (IpMod(i)-IpMod(i-1))/(IpMod(i)+IpMod(i-1)); % del_ip/2ip
    del_is(i,1) = (IsMod(i)-IsMod(i-1))/(IsMod(i)+IsMod(i-1)); % del_is/2is
    del_ipis(i,1) = (IsMod(i)+IsMod(i-1))/(IpMod(i)+IpMod(i-1)); % is/ip    
end

del_ip(isnan(del_ip)) = 0;
del_is(isnan(del_is)) = 0;
del_ipis(isnan(del_ipis)) = 0;

r = (1+tan(angle).^2).*del_ip-8.*del_ipis.^2.*sin(angle).^2.*del_is;

synthetic = conv2(r,wavelet,'same');

end

