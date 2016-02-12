function [ Vp1_r Vs1_r Rho1_r Vp2_r Vs2_r Rho2_r Por1_r Por2_r] = RockPropDist(NI,Dist,MaxZ,Pves_filt_1,Pves_filt_2,Vp1_pred,Vs1_pred,Rho1_pred,Vp2_pred,Vs2_pred,Rho2_pred,Vp1,Vs1,Rho1,Vp2,Vs2,Rho2,VpT1,VsT1,RhoT1,VpT2,VsT2,RhoT2,Por1,Por1_pred,PorT1,Por2,Por2_pred,PorT2)
%RockPropDist Summary of this function goes here
%   Calculates residuals between filtered data points and trend fits
% Uses the residuals to determine mean and covariances in Vp, Vs and Rho

Pvs=1:MaxZ;
Pves2=round(Pves_filt_2);
Pves1=round(Pves_filt_1);


for x = 1:length(Pves2);
    for y = 1:MaxZ       
        if Pvs(y) == Pves2(x)
        Res2(x,1) = (Vp2(x)-VpT2(y))./Vp2_pred;
        Res2(x,2) = (Vs2(x)-VsT2(y))./Vs2_pred;
        Res2(x,3) = (Rho2(x)-RhoT2(y))./Rho2_pred;
        Res2(x,4) = (Por2(x)-PorT2(y))./Por2_pred;
        end
    end
end

for xx = 1:length(Pves1);
    for yy = 1:MaxZ       
        if Pvs(yy) == Pves1(xx)
        Res1(xx,1) = (Vp1(xx)-VpT1(yy))./Vp1_pred;
        Res1(xx,2) = (Vs1(xx)-VsT1(yy))./Vs1_pred;
        Res1(xx,3) = (Rho1(xx)-RhoT1(yy))./Rho1_pred;
        Res1(xx,4) = (Por1(xx)-PorT1(yy))./Por1_pred;
        end
    end
end

Sigma_Res2=std(Res2);
Sigma_Res1=std(Res1);
Cov_Res2=cov(Res2);
Cov_Res1=cov(Res1); % Covariance
Mu1=[0 0 0 0];
Mu2=[0 0 0 0];

% Pick properties randomly
r1=Dist*mvnrnd(Mu1,Cov_Res1,NI);
r2=Dist*mvnrnd(Mu2,Cov_Res2,NI);

Vp1_r=Vp1_pred+(r1(:,1).*Vp1_pred);
Vs1_r=Vs1_pred+(r1(:,2).*Vs1_pred);
Rho1_r=Rho1_pred+(r1(:,3).*Rho1_pred);
Por1_r=Por1_pred+(r1(:,4).*Por1_pred);
Vp2_r=Vp2_pred+(r2(:,1).*Vp2_pred);
Vs2_r=Vs2_pred+(r2(:,2).*Vs2_pred);
Rho2_r=Rho2_pred+(r2(:,3).*Rho2_pred);
Por2_r=Por2_pred+(r2(:,4).*Por2_pred);


end
