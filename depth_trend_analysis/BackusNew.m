%Backus average%IMPORTANT:   a should not be zerofunction [Vpn,Vsn,RHOBn,Ipn,Isn,PRn]=BackusNew(Vp,Vs,RHOB,n)RHOBn=ArithmNew(RHOB,n);M=Vp.*Vp.*RHOB;G=Vs.*Vs.*RHOB;Mn=HarmonicNew(M,n);Gn=HarmonicNew(G,n);PRn=.5*(Mn./Gn-2)./(Mn./Gn-1);Ipn=sqrt(Mn.*RHOBn);Isn=sqrt(Gn.*RHOBn);Vpn=sqrt(Mn./RHOBn);Vsn=sqrt(Gn./RHOBn);