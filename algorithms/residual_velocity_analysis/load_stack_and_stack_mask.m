
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
        
trace_count = 0;
stack_all = zeros(2501,total_traces);
stack_mask_all = zeros(2501,total_traces);

for iblock = 2:23
    matfile = strcat('/segy/URY/2014_BG_water_column_imaging/matlab/rms_picks_test_block',int2str(iblock),'.mat');
    load(matfile,'stack','stack_mask','ntraces');
    trace_count = trace_count + ntraces;
    stack_all(:,trace_count-ntraces+1:trace_count) = stack;
    stack_mask_all(:,trace_count-ntraces+1:trace_count) = stack_mask;
end
    
% create vectors for scatter plot

[stky stkx] = find(stack_mask_all);
stky = stky.*2;
% decimated version

[stky10 stkx10] = find(stack_mask_all(:,10:10:end));
stky10 = stky10.*2;


