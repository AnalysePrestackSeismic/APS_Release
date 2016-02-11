function [altnodes,mcr_root,mcr_cache_root,run_script_path,uplim,lowlim] = param_node_slurm_submit
	altnodes = {'sblxpghn100'}; % {'sblxpghn100 ','sblxpgcn192 ','sblxpgcn182 '}
	mcr_root = '/apps/matlab/v2013a/'; % path location to matlab install or MCR
	mcr_cache_root = '/tmp/'; % '/localcache/mcr_cache_umask_friendly/'; % suitable temp location
	run_script_path = '/apps/gsc/matlab-library/aws_cloud_test/'; % path location to compiled matlab algorithms
	uplim=1; % lower limit on submit time
	lowlim=2; % upper limit on submit time
end
