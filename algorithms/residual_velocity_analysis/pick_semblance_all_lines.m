
% get list of lines

linedir = '/data/URY/segy/2014_BG_water_column_imaging/conditioned_datasets/semblance/';
filename = 'job_meta_';
file_ext = '.mat';

[~,filelist] = system(['find ',linedir,' -name "',filename,'*',file_ext,'" -print | sort']);

infiles = strsplit(filelist,'\n');
nlines=max(size(infiles))-1;


% for line=1:nlines;
line=45;
    linenames{line} = regexprep(infiles{line},linedir,'');
    linenames{line} = linenames{line}(1:8);
    disp([num2str(line),' ... ',linenames{line}]);
    meta_path = strcat(infiles{line});
    [alltimes alllocs allvels allfiltvels] = pick_semblance(meta_path);
    
    save(strcat(linedir,linenames{line},'/semblance_picks_',linenames{line},'.mat'));
    

% end






