
% run residuals for all blocks

% matlabpool local 8

% parfor iblock = 2:24
parfor iblock = 18:23
    
    strblock = int2str(iblock);
    iblock
    residual_velocity_analysis('/URY/segy/2014_BG_water_column_imaging/matlabjob_meta/job_meta_23Jul2014.mat',strblock,'velocity_path');
end

% matlabpool close
