%Arithmetic average (running mean) of vector 'a'
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
%with window size n
%averages n consequtive points of the vector and
%puts result to the position of the first point in the window + fix(n/2)
%adds the original vector's head and tail where the averaged values are missing

function [amean]=ArithmNew(a,n,L)

h=length(a);

for j=1:L
for i=1:h-n;
	%Mean
	Mean=sum(a(i:i+n-1,j))./n;
	amean(i+fix(n/2),j)=Mean;
end

amean(1:fix(n/2),j)=amean(fix(n/2),j);
amean(h-n+fix(n/2)+1:h,j)=a(h-n+fix(n/2)+1,j);
end

