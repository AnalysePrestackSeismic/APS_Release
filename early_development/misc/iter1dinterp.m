function [datait1di] = it1di(data, range_in, range_out)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

range_in = [((nil_stepout+1)+nil_start),(nil_end-(nil_stepout+1)),nil_shift];
range_out = [1,nil,1;1,nxl,1;1,nt,1];

[nrow,ncol,npage] = size(data);

x = (range_in(1,1):range_in(1,3):range_in(1,2));
y = (range_in(2,1):range_in(2,3):range_in(2,2));
z = (range_in(3,1):range_in(3,3):range_in(3,2));

xi = (range_out(1,1):range_out(1,3):range_out(1,2));
yi = (range_out(2,1):range_out(2,3):range_out(2,2));
zi = (range_out(3,1):range_out(3,3):range_out(3,2));

for i = 1:ncol
    for j = 1:npage
        phasei_tmp1(:,i,j) = interp1(z,data(:,i,j),zi,'linear','extrap');
    end
end
for i = 1:nt
    for j = 1:nil_block-1
        phasei_tmp2(i,:,j) = interp1(y,phasei_tmp1(i,:,j),yi,'linear','extrap');
    end
end
phasei_tmp3 = reshape(phasei_tmp2,[],(nil_block-1));
for i = 1:nt*nxl
    phasei(i,:) = interp1(z,phasei_tmp3(i,:),zi,'linear','extrap');
    if (nt*nxl/4)*floor(i/(nt*nxl/4)) == i;
        fprintf('Completed interpolation for %d / %d = %.0f%%\n',i,nt*nxl,(i*100)/(nt*nxl));
    end
end
phasei = reshape(phasei,838,651,301);
end

