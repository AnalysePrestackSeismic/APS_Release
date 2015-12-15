function [] = compile_function_extra(algorithm_name,thread_type)
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







