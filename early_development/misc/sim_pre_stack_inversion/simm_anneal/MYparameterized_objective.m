function y = MYparameterized_objective(x,Sobs,IpPRI,IsPRI,W1,W2,W3,wavelets,angles)
%function to define objective 

%extracting the Ip part
IpMOD = x(1:end/2);
%extracting the Is part
IsMOD = x(1+end/2:end);

%re-computing Smod (note that this can call out to other functions if needed)
%Smod = IpMOD*IsMOD';
strc1 = gen_smod(IpMOD,IsMOD,wavelets(:,1),angles(1));
strc2 = gen_smod(IpMOD,IsMOD,wavelets(:,2),angles(2));
strc3 = gen_smod(IpMOD,IsMOD,wavelets(:,3),angles(3));

Smod = vertcat(strc1,strc2,strc3);

% objective function
y = W1*sum(sum(abs(Sobs - Smod)))/sum(sum(abs(Sobs))) +...
        W2*(sum(abs(IpPRI - IpMOD))/sum(abs(IpPRI)) +...
            sum(abs(IsPRI - IsMOD))/sum(abs(IsPRI))) +...
                W3*( sum(abs(IsPRI./IpPRI - IsMOD./IpMOD))/sum(abs(IsPRI./IpPRI)));
