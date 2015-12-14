function [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz] = Structure2Diffusion(e1,e2,e3,l1,l2,l3,n_samp,n_xline,n_iline)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

e1 = cell2mat(e1);
e1x = reshape(e1(1,:),n_samp,n_xline,n_iline);
e1y = reshape(e1(2,:),n_samp,n_xline,n_iline);
e1z = reshape(e1(3,:),n_samp,n_xline,n_iline);

e1 = cell2mat(e2);
clear e2
e2x = reshape(e1(1,:),n_samp,n_xline,n_iline);
e2y = reshape(e1(2,:),n_samp,n_xline,n_iline);
e2z = reshape(e1(3,:),n_samp,n_xline,n_iline);

e1 = cell2mat(e3);
clear e3
e3x = reshape(e1(1,:),n_samp,n_xline,n_iline);
e3y = reshape(e1(2,:),n_samp,n_xline,n_iline);
e3z = reshape(e1(3,:),n_samp,n_xline,n_iline);
clear e1

% Construct the tensors

lambda1 = 0;
lambda2 = (2*l2.*(l2-l3))./((l1+l2).*(l2+l3));

Dxx = lambda1.*e1x.^2   + 0.5.*lambda2.*e2x.^2   + 0.5.*lambda2.*e3x.^2;
Dyy = lambda1.*e1y.^2   + 0.5.*lambda2.*e2y.^2   + 0.5.*lambda2.*e3y.^2;
Dzz = lambda1.*e1z.^2   + 0.5.*lambda2.*e2z.^2   + 0.5.*lambda2.*e3z.^2;

Dxy = lambda1.*e1x.*e1y + 0.5.*lambda2.*e2x.*e2y + 0.5.*lambda2.*e3x.*e3y;
Dxz = lambda1.*e1x.*e1z + 0.5.*lambda2.*e2x.*e2z + 0.5.*lambda2.*e3x.*e3z;
Dyz = lambda1.*e1y.*e1z + 0.5.*lambda2.*e2y.*e2z + 0.5.*lambda2.*e3y.*e3z;

end

