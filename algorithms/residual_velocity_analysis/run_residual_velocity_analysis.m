
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
% get list of filenames

datadir = '/data/URY/segy/2014_BG_water_column_imaging/conditioned_datasets/cdps/';
filescan = 'gath_block1_conditioned_cdps_s1c5_';

[infiles, nfiles] = directory_scan({datadir},filescan);

count=0;

for file = 1:nfiles
    if ~isempty(strfind(infiles.names{file},'.segy'))
        count=count+1;
        segyfiles{count}=infiles.names{file};
    end
end

nlines=size(segyfiles,2);

% for file = 1:nlines
%     infile = infiles.names{file};
%     linename{file} = regexprep(infile,'gath_block1_conditioned_cdps_s1c5_','');
%     linename{file} = regexprep(linename{file},'.segy','');
%     datadir2 = strcat(datadir,linename{file},'/');
%     segy_make_job(datadir2,infile,'189','21','37','0','0',datadir2);
% end

% run the residual velocities
% node_slurm_submit runs all blocks in each meta file
% -> first and last will fail but doesn't matter

% parfor ii = 1:nlines
for ii = 1:nlines
    infile = segyfiles{ii};
    linename{ii} = regexprep(infile,'gath_block1_conditioned_cdps_s1c5_','');
    linename{ii} = regexprep(linename{ii},'.segy','');
    meta_path = strcat(datadir,'/',linename{ii},'/job_meta/job_meta_13Oct2014.mat');
    outdir = strcat(datadir,'/',linename{ii},'/');
    gather_meta = load(meta_path);
    
    node_slurm_submit('residual_velocity_analysis',meta_path,'UK1','1',outdir,'1');
    
%     blocks = gather_meta.liveblocks(2:(end-1));
%     for iblock = blocks(1):blocks(end)
%         residual_velocity_analysis(meta_path,int2str(iblock),outdir,1);
%     end
end

