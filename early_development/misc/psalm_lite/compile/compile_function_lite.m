function [] = compile_function_lite(algorithm_name)
% function to compile the specified matlab function

func_path = strcat('/apps/gsc/matlab-mcode-beta/eslib/psalm_lite/algorithms/',algorithm_name,'/');

%     mcc_call = sprintf('mcc -v -M ''-O2'' -o %s -W main:%s -T link:exe -d %s -R ''-nojvm,-nodisplay,-singleCompThread'' %s.m',...
%      algorithm_name,algorithm_name,func_path,algorithm_name');

    mcc_call = sprintf('mcc -v -M ''-O2'' -o %s -W main:%s -T link:exe -d %s -R ''-nojvm,-nodisplay'' %s.m',...
     algorithm_name,algorithm_name,func_path,algorithm_name'); 
 
%   mcc_call = sprintf('mcc -v -o %s -W main:%s -T link:exe -d %s -R ''-nojvm,-nodisplay,'' %s.m',...
%    algorithm_name,algorithm_name,func_path,algorithm_name');

eval(mcc_call);

fprintf('Function %s successfully compiled to:\n%s\n',algorithm_name,func_path);

end







