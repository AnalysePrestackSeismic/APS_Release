

for i_file = 1:1:38
    
node_slurm_submit_lite('trim_calculation',1,140,['/data/KEN/segy/2012_L10ab_final/raw_migration_gathers/109j05_post-mig_gather_segy_disk_sw',num2str(i_file),'.mat_lite'],'UK1','nodelistcj','2','segy','/data/KEN/segy/2012_L10ab_final/james_trim_gathers/','4','5','25','ad');   
% srun -p UK1 -c 10 -J trim_calculation_block_1 /apps/gsc/matlab-mcode-beta/eslib/psalm_lite/algorithms/run_function.sh /apps/matlab/v2011a trim_calculation /data/KEN/segy/2012_L10ab_final/raw_migration_gathers/109j05_post-mig_gather_segy_disk_sw15.mat_lite 1 110 segy /data/KEN/segy/2012_L10ab_final/james_trim_gathers/ 4 5 25 ad &    
%     for i_block = 1:1:250
%         trim_calculation(['/data/KEN/segy/2012_L10ab_final/raw_migration_gathers/109j05_post-mig_gather_segy_disk_sw',num2str(i_file),'.mat_lite'],num2str(i_block),'250','','/data/NOR/segy/hegg_fwi_test_2012/gathers/output/','4','5','25','{}')
%     end
end