function [] = add_including_poly_to_job_meta(job_meta_path,ils,xls)

job_meta.ply.il=ils;
job_meta.poly.xl=xls;
poly.method=1;
ils=double(ils);
xls=double(xls);
ils=ils';
xls=xls';


inline_sort_blocks(job_meta_path)
job_meta = load(job_meta_path);                 % Load job meta information
counter=1;

for i_block=1:length(job_meta.liveblocks);
    blk=job_meta.liveblocks(i_block);
    ilxl_blk=job_meta.block_keys(blk,:);
    il_corner_pts= [ilxl_blk(1) ilxl_blk(2)];
    xl_corner_pts= [ilxl_blk(3) ilxl_blk(4)];
    
    [IN ON]= inpolygon(il_corner_pts,xl_corner_pts,ils,xls);
    if sum(IN)==2
       liveblocks_poly(counter)=job_meta.liveblocks(i_block);
       counter=counter+1;
    end
end
job_meta.liveblocks = liveblocks_poly';
str_date = date;
str_date = regexprep(str_date, '-', '');
job_meta_path_poly=strcat(job_meta.output_dir,'job_meta/','job_meta_poly_',str_date,'.mat');

save(job_meta_path_poly,'-struct','job_meta','-v7.3');
end
