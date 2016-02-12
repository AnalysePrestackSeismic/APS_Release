%function [filename,filename_out,md_rs,age_rs] = LogInterpolation(well,md,age)
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
%LOGINTERPOLATION Summary of this function goes here
%   Detailed explanation goes here

filename='Chaza1_md_age.txt';
filename_out='Chaza1_md_age_interpolated.txt';
[well md age]=textread(filename,'%s %f %f');
type='linear';

min_val=min(md);
max_val=max(md);

delta=0.1524;

md_rs=interp1(md,md,min_val:delta:max_val,type)';
age_rs=interp1(md,age,min_val:delta:max_val,type)';
out=[md_rs age_rs];
dlmwrite(filename_out,out,'delimiter',' ','precision',8);

%end

