function [] = crop_live_blocks(job_meta_path,il1,xl1,il2,xl2)

inline_sort_blocks(job_meta_path)
job_meta = load(job_meta_path);                 % Load job meta information
counter=1;

for i_block=1:length(job_meta.liveblocks);
    blk=job_meta.liveblocks(i_block);
    ilxl_blk=job_meta.block_keys(blk,:);
    if (ilxl_blk(1)>il1 && ilxl_blk(2)<il2 && ilxl_blk(3)>xl1 && ilxl_blk(4)<xl2)
       liveblocks_crop(counter)=job_meta.liveblocks(i_block);
       counter=counter+1;
    end
end
job_meta.liveblocks = liveblocks_crop';
str_date = date;
str_date = regexprep(str_date, '-', '');
job_meta_path_crop=strcat(job_meta.output_dir,'job_meta/','job_meta_crop',str_date,'.mat');

save(job_meta_path_crop,'-struct','job_meta','-v7.3');
end
