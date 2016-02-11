function segy_plot_live_blocks_aperture(job_meta_path,i_vol,aperture)
%% ------------------ Disclaimer  ------------------
% 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) makes no representation or warranty, express or implied, in 
% respect to the quality, accuracy or usefulness of this repository. The code
% is this repository is supplied with the explicit understanding and 
% agreement of recipient that any action taken or expenditure made by 
% recipient based on its examination, evaluation, interpretation or use is 
% at its own risk and responsibility.
% 
% No representation or warranty, express or implied, is or will be made in 
% relation to the accuracy or completeness of the information in this 
% repository and no responsibility or liability is or will be accepted by 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) in relation to it.
%% ------------------ License  ------------------ 
% GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
%% github
% https://github.com/AnalysePrestackSeismic/
%% ------------------ FUNCTION DEFINITION ---------------------------------
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