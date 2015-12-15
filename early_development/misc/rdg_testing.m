function [] = rdg_testing(job_meta_path)

job_meta=load(job_meta_path);

for k=1:520
    bl=job_meta.liveblocks(k);
    bls=num2str(bl);
    load_wavelet_spatial('/data/KEN/segy/2013_L10AB/depth/final_stacks/from_usb/output_depth/job_meta/job_meta_12Jun2014.mat', '/data/KEN/segy/2013_L10AB/depth/final_stacks/from_usb/output_depth/wavelets/',bls,'2');
    fprintf('Passed: %g \n',bl);
    
end

end