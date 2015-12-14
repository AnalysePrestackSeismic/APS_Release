function segy_plot_blocks(job_meta_path)

job_meta = load(job_meta_path);

figure

% For this example, there are five (m=5) triangles (n=3). 
% The total number of vertices is mxn, or k = 15. 
% xdata = [2 2 0 2 5;
%          2 8 2 4 5;
%          8 8 2 4 8];
% ydata = [4 4 4 2 0;
%          8 4 6 2 2;
%          4 0 4 0 0];
% zdata = ones(3,5);
% p = patch(xdata,ydata,zdata,'b')
% 
% clear cdata 
% cdata(:,:,1) = [0 0 1 0 0.8];
% cdata(:,:,2) = [0 0 0 0 0.8];
% cdata(:,:,3) = [1 1 1 0 0.8];
% set(p,'FaceColor','flat','CData',cdata)




% should we initialize the variables?

for i_block = 1:1:size(job_meta.block_keys,1)
cjxdata(:,i_block) = [job_meta.block_keys(i_block,1); job_meta.block_keys(i_block,2); job_meta.block_keys(i_block,2); job_meta.block_keys(i_block,1)];
cjydata(:,i_block) = [job_meta.block_keys(i_block,3); job_meta.block_keys(i_block,3); job_meta.block_keys(i_block,4); job_meta.block_keys(i_block,4)];

cdata(1,i_block,1) = 1;
cdata(1,i_block,2) = 1;
cdata(1,i_block,3) = 1;

end
%zdata = ones(4,size(job_meta.block_keys,1));
%patch(cjxdata,cjydata,zdata,'w');
p = patch(cjxdata,cjydata,'w');
%set(p,'FaceColor','flat','CData',cdata)

%hold all;
%set(p,'FaceColor',[0 0.1 0]);
loopfin = size(job_meta.liveblocks,1);
lpi = 1;
while lpi <= loopfin
i_block = job_meta.liveblocks(lpi);
%cjbxdata(:,lpi) = [job_meta.block_keys(i_block,1); job_meta.block_keys(i_block,2); job_meta.block_keys(i_block,2); job_meta.block_keys(i_block,1)];
%cjbydata(:,lpi) = [job_meta.block_keys(i_block,3); job_meta.block_keys(i_block,3); job_meta.block_keys(i_block,4); job_meta.block_keys(i_block,4)];

cdata(1,i_block,1) = 0;
cdata(1,i_block,2) = 1;
cdata(1,i_block,3) = 0;

lpi = lpi + 1;
end
%
set(p,'FaceColor','flat','CData',cdata)

%patch(cjbxdata,cjbydata,'g');

%patch(cjbxdata,cjbydata,'w','Facecolor',[0,0.1,0]);

% ascii read for uncompleted blocks
%stuckjob = dlmread('/home/jaffa/jonesce/stuckjobs');

%loopfin = size(job_meta.liveblocks,1);
%lpi = 1;

% 
% [seismic, ~, ~, ~] = node_segy_read(job_meta_path,i_vol,'1');
% 
% figure
% 
% scatter(seismic.trace_ilxl_bytes(:,1),seismic.trace_ilxl_bytes(:,2))
% hold all
% scatter(seismic.trace_ilxl_bytes(:,1),seismic.trace_ilxl_bytes(:,4)) 



% 
% scatter([job_meta.block_keys(1,1);job_meta.block_keys(end,1)],[job_meta.block_keys(1,3);job_meta.block_keys(end,3)])
% 
% for i_block = 1:1:size(job_meta.block_keys,1)
%     colour = [0.98,0.78,0.49];
%     colour = [random('bino',1,0.5),random('bino',1,0.5),random('bino',1,0.5)];
%     fill([job_meta.block_keys(i_block,2);job_meta.block_keys(i_block,2);...
%         job_meta.block_keys(i_block,1);job_meta.block_keys(i_block,1)],...
%         [job_meta.block_keys(i_block,3);job_meta.block_keys(i_block,4);...
%         job_meta.block_keys(i_block,4);job_meta.block_keys(i_block,3)],colour)
% end

%  while lpi <= loopfin
%      i_block = job_meta.liveblocks(lpi);
%      %colour = [random('bino',1,0.5),random('bino',1,0.5),random('bino',1,0.5)];
%      %colour = [0.77,0.96,0.73];
%      colour = [0.98,0.78,0.49];
%      fill([job_meta.block_keys(i_block,2);job_meta.block_keys(i_block,2);...
%          job_meta.block_keys(i_block,1);job_meta.block_keys(i_block,1)],...
%          [job_meta.block_keys(i_block,3);job_meta.block_keys(i_block,4);...
%          job_meta.block_keys(i_block,4);job_meta.block_keys(i_block,3)],colour)
%      lpi = lpi + 1;
%  end
% 
%  
% scatter(seismic.trace_ilxl_bytes(:,1),seismic.trace_ilxl_bytes(:,2))
% hold all
% scatter(seismic.trace_ilxl_bytes(:,1),seismic.trace_ilxl_bytes(:,4)) 
 
% while lpi <= size(stuckjob,1)
%     i_block = stuckjob(lpi);
%     %colour = [random('bino',1,0.5),random('bino',1,0.5),random('bino',1,0.5)];
%     %colour = [0.77,0.96,0.73];
%     colour = [0.9,0.08,0.08];
%     fill([job_meta.block_keys(i_block,2);job_meta.block_keys(i_block,2);...
%         job_meta.block_keys(i_block,1);job_meta.block_keys(i_block,1)],...
%         [job_meta.block_keys(i_block,3);job_meta.block_keys(i_block,4);...
%         job_meta.block_keys(i_block,4);job_meta.block_keys(i_block,3)],colour)
%     lpi = lpi + 1;
% end



end