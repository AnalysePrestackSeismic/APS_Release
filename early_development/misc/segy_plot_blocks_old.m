function segy_plot_blocks_old()




runtime = java.lang.Runtime.getRuntime();
process = runtime.exec('program arg1 arg2');  % non-blocking
% Continue Matlab processing in parallel to spawned process
%When we need to collect scalar results, we could use the processâ€™ result code:
rc = process.waitFor();    % block Matlab until external program ends
rc = process.exitValue();  % fetch an ended process' return code
%Or, if we need to abandon the work, we could stop the spawned process:
process.destroy();         % force-kill the process (rc will be 1)













end