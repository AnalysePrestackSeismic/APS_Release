function [] = compile_function_extra(algorithm_name,thread_type)
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
% COMPILE_FUNCTION: to compile the specified matlab function
% Inputs:
%  - algorithm_name - self explanatory. Algorithms should be stored in a
%  subdirectory of the algorithms folder with the same name as the function
%  to be compiled.
%  - thread_type - 1 = Compile single thread
%                  2 = Compile multi threaded

thread_type = str2double(thread_type);

func_path = strcat('/apps/gsc/matlab-library/development/maps/extra/');

if thread_type == 1
    % Single thread
    mcc_call = sprintf('mcc -v -M ''-O2'' -o %s -W main:%s -T link:exe -d %s -R ''-nojvm,-nodisplay,-singleCompThread'' %s.m',...
        algorithm_name,algorithm_name,func_path,algorithm_name');
    
elseif thread_type == 2
    % Multithread
    mcc_call = sprintf('mcc -v -o %s -W main:%s -T link:exe -d %s -R ''-nojvm,-nodisplay,'' %s.m',...
        algorithm_name,algorithm_name,func_path,algorithm_name');
end

eval(mcc_call);

fprintf('Function %s successfully compiled to:\n%s\n',algorithm_name,func_path);

end







