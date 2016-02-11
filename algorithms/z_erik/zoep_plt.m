function zoep_plt(rho1,a1,b1,rho2,a2,b2,incwav,anginc,flag)
%% ------------------ Disclaimer  ------------------
% 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) makes no representation or warranty, express or implied, in 
% respect to the quality, accuracy or usefulness of this repository. The code
% is this repository is supplied with the explicit understanding and 
% agreement of recipient that any action taken or expenditure made by 
% recipient based on its examination, evaluation, interpretation or use is 
% at its own risk and responsibility.
% 
% No representation or warranty, express or implied, is or will be made in 
% relation to the accuracy or completeness of the information in this 
% repository and no responsibility or liability is or will be accepted by 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) in relation to it.
%% ------------------ License  ------------------ 
% GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
%% github
% https://github.com/AnalysePrestackSeismic/
%% ------------------ FUNCTION DEFINITION ---------------------------------
%
% zoepplt(rho1,a1,b1,rho2,a2,b2,incwav,anginc,flag)
% zoepplt(rho1,a1,b1,rho2,a2,b2,incwav,anginc)
% zoepplt(rho1,a1,b1,rho2,a2,b2,incwav)
%
% ZOEPPLT computes displacement reflection and transmission coefficients
% using the ZOEPPRITZ function. Zoepplt produces a plot of these 
% coefficients vs. angle of incidence, and will show all critical
% angles as dashed lines if so desired. 
%
% rho1    = density of the incidence medium
% a1      = p-wave velocity of incidence medium
% b1      = s-wave velocity of incidence medium
% rho2    = density of the transmission medium
% a2      = p-wave velocity of transmission medium
% b2      = s-wave velocity of transmission medium
% incwav  = 1 for incident p-wave 
%	  = 2 for incident s-wave
% anginc  = incidence angles (in degrees)
% ************ Default = 0 to 90 in steps of 1 degree *******************
%
% flag   = 1 shows the critical angle(s) as a dashed line on plot
%        = 2 shows critical angles plus numeric value
%	 = 3 does not show critical angle(s)
%************* Default = 1 (shows critical angles) **********************
%  
% NOTE: if b1 and b2 = 0, a liquid-liquid interface is treated
% NOTE: if b1 = 0 and b2 is non-zero, a liquid-solid interface is treated 
% NOTE: if b1 is non-zero and b2 = 0, a solid-liquid interface is treated
% NOTE: if a1 and a2 = 0, the incident SH wave case is treated
% 
% by Jeff Larsen, CREWES Project, May 1997
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

RCP='P wave reflection coefficient vs. angle of incidence';
RCS='S wave reflection coefficient vs. angle of incidence';
TCP='P wave transmission coefficient vs. angle of incidence';
TCS='S wave transmission coefficient vs. angle of incidence';
TC='TC';
RC='RC';
ANGINC='Angle of Incidence (degrees)';
inc=1;

if nargin<8
  anginc=0:90;
end

if nargin<9
  flag=2;
end

hfig = figcent(.5,.9);  
set(hfig,'menubar','none','units','normalized');
hmenu = uimenu(gcf,'label','Actions');
  hzoom = uimenu(hmenu,'label','Zoom','callback','simplezoom');

%================= Find the critical angles ==============================
if flag==1 | flag==2
	dpr = 180/pi;
	ca1=0; ca2=0; cb2=0; critat=0; critbt=0; critar=0;

	if incwav==1
		velincmed = a1;
	elseif incwav==2
		velincmed = b1;
	end
	if a1~=0
		ca1 = (asin(velincmed/a1))*dpr;
	end
	if a2~=0
		ca2 = (asin(velincmed/a2))*dpr;
	end
	if b2~=0
		cb2 = (asin(velincmed/b2))*dpr; 
	end
	if (real(ca1)==90)
		ca1 = 0;
	end
	if (real(ca2)==90)
		ca2 = 0;
	end
	if (real(cb2)==90)
		cb2 = 0;
	end
else 
	ca1=0; ca2=0; cb2=0; 
end
%============ SH Wave case for a Solid-Solid Interface ==============
if (a1==0)&(a2==0)&(b1>0)&(b2>0)
  if incwav==2
 	newtitle='SH wave case for a Solid-Solid Interface';
  
 	for irfwav=2:2:4
  	coef=zoeppritz(rho1,a1,b1,rho2,a2,b2,incwav,irfwav,1,anginc);
  	coef=real(coef);
  	subplot(2,1,inc), plot(anginc,coef);
  	grid;
  	xlabel(ANGINC);
  	axis([min(anginc) max(anginc) min(coef) max(coef)]);
    		if irfwav==2
     		 title(RCS,'VerticalAlignment','cap','fontweight','bold');
     		 ylabel(RC);
      		elseif irfwav==4
      		title(TCS,'VerticalAlignment','cap','fontweight','bold');
      		ylabel(TC);
      		if( cb2 )
      		   line([cb2 cb2],[min(coef) max(coef)],'color','r','LineStyle','--');
      		   if( flag==2 )
             	     text(cb2,min(coef),num2str(cb2),'color',[0.6 0 .6], ...
                     'HorizontalAlignment','center','fontweight','bold');
                   end
      		end
    		end
  	inc=inc+1;
  	end
  elseif incwav==1
  	close;
  	error('You must have an incident S wave'); 
  end
%=================== Solid-Solid Interface ==========================
elseif (b1>0)&(b2>0)&(a1>0)&(a2>0)
  
  if incwav==1
  	newtitle='Solid-Solid Interface for an incident P wave';
  elseif incwav==2
  	newtitle='Solid-Solid Interface for an incident S wave';
  end
    
  for irfwav=1:4
    coef=zoeppritz(rho1,a1,b1,rho2,a2,b2,incwav,irfwav,1,anginc);
    coef=real(coef);
    
    subplot(4,1,irfwav), plot(anginc,coef);
      if irfwav==1
      	   if( ca1 )
     	      line([ca1 ca1],[min(coef) max(coef)],'color','r','LineStyle','--');
     	      if( flag==2 )
       	        text(ca1,min(coef),num2str(ca1),'color',[0.6 0 .6], ...
                'HorizontalAlignment','center','fontweight','bold');
              end
     	   end
     	  title(RCP,'VerticalAlignment','cap','fontweight','bold');
      	  ylabel(RC);
      elseif irfwav==2
      	  title(RCS,'VerticalAlignment','cap','fontweight','bold'); 
      	  ylabel(RC);
      elseif irfwav==3
         if( ca2 )
       	    line([ca2 ca2],[min(coef) max(coef)],'color','r','LineStyle','--');
       	    if( flag==2 )
       	      text(ca2,min(coef),num2str(ca2),'color',[0.6 0 .6], ...
              'HorizontalAlignment','center','fontweight','bold');
            end
       	 end
         title(TCP,'VerticalAlignment','cap','fontweight','bold');
         ylabel(TC);
      elseif irfwav==4
         if( cb2 )
            line([cb2 cb2],[min(coef) max(coef)],'color','r','LineStyle','--');
            if( flag==2 )
              text(cb2,min(coef),num2str(cb2),'color',[0.6 0 .6], ...
              'HorizontalAlignment','center','fontweight','bold');
            end
         end
         title(TCS,'VerticalAlignment','cap','fontweight','bold');
         ylabel(TC);
      end
    grid;
    xlabel(ANGINC);
    axis([min(anginc) max(anginc) min(coef) max(coef)]);
    
    end
%====================================================================
orient tall; % gives better printing result
else
	close;
	error('Please check your input parameters again');	
end
set(hfig,'name',newtitle,'units','pixels');
	
 	
