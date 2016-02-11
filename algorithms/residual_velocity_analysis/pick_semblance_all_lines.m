
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






