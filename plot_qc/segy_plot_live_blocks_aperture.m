function segy_plot_live_blocks_aperture(job_meta_path,i_vol,aperture)

job_meta = load(job_meta_path);

job_meta.block_keys(:,1) = job_meta.block_keys(:,1)-str2num(aperture);
job_meta.block_keys(:,2) = job_meta.block_keys(:,2)+str2num(aperture);
job_meta.block_keys(:,3) = job_meta.block_keys(:,3)-str2num(aperture);
job_meta.block_keys(:,4) = job_meta.block_keys(:,4)+str2num(aperture);

figure


for i_block = 1:1:size(job_meta.block_keys,1)
cjxdata(:,i_block) = [job_meta.block_keys(i_block,1); job_meta.block_keys(i_block,2); job_meta.block_keys(i_block,2); job_meta.block_keys(i_block,1)];
cjydata(:,i_block) = [job_meta.block_keys(i_block,3); job_meta.block_keys(i_block,3); job_meta.block_keys(i_block,4); job_meta.block_keys(i_block,4)];

cdata(1,i_block,1) = 1;
cdata(1,i_block,2) = 1;
cdata(1,i_block,3) = 1;

end
%zdata = ones(4,size(job_meta.block_keys_aperture,1));
%patch(cjxdata,cjydata,zdata,'w');
p = patch(cjxdata,cjydata,'w');
%set(p,'FaceColor','flat','CData',cdata)

%hold all;
%set(p,'FaceColor',[0 0.1 0]);
loopfin = size(job_meta.liveblocks,1);
lpi = 1;
while lpi <= loopfin
i_block = job_meta.liveblocks(lpi);
%cjbxdata(:,lpi) = [job_meta.block_keys_aperture(i_block,1); job_meta.block_keys_aperture(i_block,2); job_meta.block_keys_aperture(i_block,2); job_meta.block_keys_aperture(i_block,1)];
%cjbydata(:,lpi) = [job_meta.block_keys_aperture(i_block,3); job_meta.block_keys_aperture(i_block,3); job_meta.block_keys_aperture(i_block,4); job_meta.block_keys_aperture(i_block,4)];

cdata(1,i_block,1) = 0;
cdata(1,i_block,2) = 0.95;
cdata(1,i_block,3) = 0;

lpi = lpi + 1;
end
%
set(p,'FaceColor','flat','CData',cdata)

end