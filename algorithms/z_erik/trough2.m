inter=[1,2,3,3,2,1,1,2,3,3,2,1,1,2,3]
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
deriv_1_inter=diff(inter)
inter_t=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
grad_t=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
sign_d1=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

for i=2:14
    if deriv_1_inter(i)>0
        sign_d1(i)=1;
    elseif deriv_1_inter(i)<0
        sign_d1(i)=-1;
    end
    if sign_d1(i-1)<sign_d1(i)
            inter_t(i)=inter(i);
			grad_t(i)=grad(i);
    end
end

sign_d1
inter_t
grad_t