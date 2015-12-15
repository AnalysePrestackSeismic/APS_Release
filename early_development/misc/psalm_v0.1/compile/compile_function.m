function [] = compile_function(algorithm_name)
% function to compile the specified matlab function

func_path = strcat('/apps/gsc/matlab-library/psalm_v0.1/algorithms/',algorithm_name);

mcc_call = sprintf('mcc -o %s -W main:%s -T link:exe -d %s -R ''-nojvm,-nodisplay,-singleCompThread'' %s.m',...
    algorithm_name,algorithm_name,func_path,algorithm_name');

eval(mcc_call);

fprintf('Function %s successfully compiled to:\n%s\n',algorithm_name,func_path);

% mcc -o test -W main:test -T link:exe -d /home/thedailyrant/compiler/test/src 
% -w enable:specified_file_mismatch -w enable:repeated_file 
% -w enable:switch_ignored -w enable:missing_lib_sentinel 
% -w enable:demo_license -f /usr/bin/gcc -R -nojvm -R -nodisplay 
% -R '-R -singleCompThread' % 
% -v /home/thedailyrant/mcode/algorithms/ig_inversion/synang.m 
% -a /home/thedailyrant/mcode/node/node_make_processing_positions.m 
% -a /home/thedailyrant/mcode/node/node_segy_read_traces.m 
% distribute for slice or trace

end