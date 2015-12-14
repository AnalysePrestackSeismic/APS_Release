function [t,r]=shootrayvxz_bg(tstep,r0)
% SHOOTRAYVXZ RK4 raytracing in v(x,z) with nearest neighbor int.
%
% [t,r]=shootrayvxz(tstep,r0)
%
% Ray tracer for a 2D gridded velocity model. The model
% is built by defining a matrix of velocity values and calling
% RAYVELMOD. This establishes 3 global matrices containing v^2
% and the logarithmic derivatives of v with respect to x and z.
% (The velocity model matrix may be deleted after calling RAYVELMOD
% if memory is limited.) The raytracer implements an RK4 (4th order 
% Runge-Kutta) solution to the ray equations and, by default, uses
% nearest neighbor interpolation. Bilinear interpolation is available
% for more accuracy. To get bilinear interpolation, use SHOOTRAYVXZ_G. 
% Rays will be terminated when they reach the maximum time or leave the 
% bounds of the velocity model.
%
% tstep ... vector of time steps, (e.g [0:.01:tmax] steps from 0 to
%           tmax in 10 millisecond intervals.)
% r0 ... initial values of the ray vector. r0(1) and r0(2) are the
%        starting x and z coordinates and r0(3) and r0(4) are the
%        horizontal and vertical slownesses.
% t ... output time vector. If the ray stays within the bounds of
%        the model, this is the same as tstep, otherwise it may be
%        smaller.
% r ... output ray matrix. This is an N by 4 matrix, where N=length(t),
%       with each row being a ray vector for the corresponding time. The
%       ray vectors are as described in r0 above.
%
% G.F. Margrave and P.F. Daley, CREWES, June 2000
%
% NOTE: It is illegal for you to use this software for a purpose other
% than non-profit education or research UNLESS you are employed by a CREWES
% Project sponsor. By using this software, you are agreeing to the terms
% detailed in this software's Matlab source file.
 
% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by 
% its author (identified above) and the CREWES Project.  The CREWES 
% project may be contacted via email at:  crewesinfo@crewes.org
% 
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) Use of this SOFTWARE by any for-profit commercial organization is
%    expressly forbidden unless said organization is a CREWES Project
%    Sponsor.
%
% 2) A CREWES Project sponsor may use this SOFTWARE under the terms of the 
%    CREWES Project Sponsorship agreement.
%
% 3) A student or employee of a non-profit educational institution may 
%    use this SOFTWARE subject to the following terms and conditions:
%    - this SOFTWARE is for teaching or research purposes only.
%    - this SOFTWARE may be distributed to other students or researchers 
%      provided that these license terms are included.
%    - reselling the SOFTWARE, or including it or any portion of it, in any
%      software that will be resold is expressly forbidden.
%    - transfering the SOFTWARE in any form to a commercial firm or any 
%      other for-profit organization is expressly forbidden.
%
% END TERMS OF USE LICENSE

global RTV2 RTDLNVDX RTDLNVDZ RTDG RTDX RTX RTZ
% global BILINEAR

if(isempty(RTDG))
	error('velocity model not defined. Use RAYVELMOD')
end


[m,n]=size(RTV2);                   % Size of velocity model, remember the velocity model is padded in 
edgex=4;                            % Abort within this many dg of edge
topx=2;
botx=4;
xmax=RTX(n-edgex);
zmax=RTZ(m-botx);
xmin=RTX(edgex);
zmin=RTZ(topx);

r0=r0(:)';                          % The initial ray vector

r=zeros(length(tstep),length(r0));  % Initialize the matricx to store initial ray vector
r(1,:)=r0;                          % The first point in the ray matrix is the initial ray vector supplied

for k=2:length(tstep)
% 	tnow=tstep(k-1);
	rnow=r(k-1,:);
	h=(tstep(k)-tstep(k-1));

   %k1
   xx=(rnow(1)-RTX(1))/RTDX+1; zz=(rnow(2)-RTZ(1))/RTDG+1;% Actual fractional grid
   
   if round(xx)>n || round(zz)>m
       r(k:length(tstep),:)=repmat(rnow,(length(tstep)-k+1),1);
       break;
   end
   
   %nearest neighbor interpolation
   v2=RTV2(round(zz),round(xx));
   dlnvdx=RTDLNVDX(round(zz),round(xx));
   dlnvdz=RTDLNVDZ(round(zz),round(xx));
   k1=h*[v2*rnow(3),v2*rnow(4),-dlnvdx,-dlnvdz];

   %k2
   rk=rnow+.5*k1;
   xx=(rk(1)-RTX(1))/RTDX+1;zz=(rk(2)-RTZ(1))/RTDG+1;
   if round(xx)>n || round(zz)>m
       r(k:length(tstep),:)=repmat(rnow,(length(tstep)-k+1),1);
       break;
   end
   v2=RTV2(round(zz),round(xx));
   dlnvdx=RTDLNVDX(round(zz),round(xx));
   dlnvdz=RTDLNVDZ(round(zz),round(xx));
   k2=h*[v2*rk(3),v2*rk(4),-dlnvdx,-dlnvdz];

   %k3
   rk=rnow+.5*k2;
   xx=(rk(1)-RTX(1))/RTDX+1;zz=(rk(2)-RTZ(1))/RTDG+1;
   if round(xx)>n || round(zz)>m
       r(k:length(tstep),:)=repmat(rnow,(length(tstep)-k+1),1);
       break;
   end
   v2=RTV2(round(zz),round(xx));
   dlnvdx=RTDLNVDX(round(zz),round(xx));
   dlnvdz=RTDLNVDZ(round(zz),round(xx));
   k3=h*[v2*rk(3),v2*rk(4),-dlnvdx,-dlnvdz];

   %k4
   rk=rnow+k3;
   xx=(rk(1)-RTX(1))/RTDX+1;zz=(rk(2)-RTZ(1))/RTDG+1;
   if round(xx)>n || round(zz)>m
       r(k:length(tstep),:)=repmat(rnow,(length(tstep)-k+1),1);
       
       break;
   end
   v2=RTV2(round(zz),round(xx));
   dlnvdx=RTDLNVDX(round(zz),round(xx));
   dlnvdz=RTDLNVDZ(round(zz),round(xx));
   k4=h*[v2*rk(3),v2*rk(4),-dlnvdx,-dlnvdz];

   %solve
	r(k,:)=rnow+k1/6+k2/3+k3/3+k4/6;
	
	%renormalize p and q
	vapp=(r(k,3)^2+r(k,4)^2)^(-.5);
	r(k,3)=r(k,3)*vapp/sqrt(v2);
	r(k,4)=r(k,4)*vapp/sqrt(v2);

 	if((r(k,1)>xmax && r(k,3)>0) || (r(k,1)<xmin && r(k,3)<0) || r(k,2)>zmax || r(k,2)<zmin)
 		break
 	end

end

if(k<length(tstep))
   r(k+1:length(tstep),:)=[];
   t=tstep(1:k);
else
   t=tstep;
end
