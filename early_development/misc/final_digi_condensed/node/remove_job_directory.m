function remove_job_directory(job_dir)

job_dir = sprintf('rm -rf %s%s',job_dir);
system(job_dir);

end