load('/apps/gsc/matlab-library/development/matlab_comp_test.mat');
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
dsign_1_deriv = sign(dmax_1st_deriv);
dsign_1_deriv(dsign_1_deriv == 0) = 1;
ddiffsign = diff(dsign_1_deriv);
dmax_1st_deriv = diff(newposvals);
dsign_1_deriv = sign(dmax_1st_deriv);
dsign_1_deriv(dsign_1_deriv == 0) = 1;
ddiffsign = diff(dsign_1_deriv);
figure(1); plot(origuprloc,newposvals,'-xb'); hold on ; plot(origuprloc(uplocs),newposvals(uplocs),'om');  plot(origuprloc(downlocs),newposvals(downlocs),'og'); hold off;


alllocs = seqno([-1;ddiffsign;-1] < 0 | [1;ddiffsign;1] > 0 );

figure(2); plot(origuprloc,newposvals,'-xb'); hold on ; plot(origuprloc(alllocs),newposvals(alllocs),'om');

newinterp = double(makefastinterp1(double(origuprloc(alllocs)),double(newposvals(alllocs)),double(origuprloc)));


figure(3); plot(origuprloc,newposvals,'-xb'); hold on ; plot(origuprloc(alllocs),newposvals(alllocs),'om');plot(origuprloc,newinterp,'-og');

diffdiff = cumsum(double(newinterp - newposvals));
difftoadd = [1;diff(diffdiff(downlocs),1)]';


figure(3); plot(origuprloc,newposvals,'-xb'); hold on ; plot(origuprloc(alllocs),newposvals(alllocs),'om');plot(origuprloc,newinterp,'-og');plot(origuprloc(downlocs),difftoadd,'or'); hold off


 figure(3); plot(origuprloc,newposvals,'-xb'); hold on ; plot(origuprloc(alllocs),newposvals(alllocs),'om');plot(origuprloc,newinterp,'-og');plot(origuprloc(uplocs),difftoadd,'or'); hold off


 
 figure(3); plot(origuprloc,newposvals,'-xb'); hold on ; plot(origuprloc(alllocs),newposvals(alllocs),'xm');plot(origuprloc,newinterp,'-og');plot(origuprloc(uplocs),newposvals(uplocs)-difftoadd','or'); hold off



difftoadd = [1;diff(diffdiff(downlocs),1)];
difftoadd(difftoadd<0) = 0;
newposvalsold = newposvals;
newposvals(uplocs) =  newposvals(uplocs)-difftoadd;