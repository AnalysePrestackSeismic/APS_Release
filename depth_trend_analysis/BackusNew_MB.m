%Backus average%IMPORTANT:   a should not be zerofunction [Vpn,Vsn,RHOBn]=BackusNew_MB(Vp,Vs,RHOB,n)RHOBn=ArithmNew_MB(RHOB,n,1);M=Vp.*Vp.*RHOB;G=Vs.*Vs.*RHOB;Mn=HarmonicNew_MB(M,n);Gn=HarmonicNew_MB(G,n);PRn=.5*(Mn./Gn-2)./(Mn./Gn-1);Ipn=sqrt(Mn.*RHOBn);Isn=sqrt(Gn.*RHOBn);Vpn=sqrt(Mn./RHOBn);Vsn=sqrt(Gn./RHOBn);