
% scan sailline datasets

% get list of filenames

datadir = '/data/URY/segy/2014_BG_water_column_imaging/conditioned_datasets/semblance/';
filescan = 'gath_block1_semblance_io_s1c5_';

[infiles, nfiles] = directory_scan({datadir},filescan);

for file = 1:nfiles
    infile = infiles.names{file};
    linename = regexprep(infile,'gath_block1_semblance_io_s1c5_','');
    linename = regexprep(linename,'.segy','');
    datadir2 = strcat(datadir,linename,'/');
    segy_make_job(datadir2,infile,'189','193','37','0','0',datadir2);
end
