% Search for best model/checkshot match
%
% Need to look for most consistent error in some way
%

%% Set up variables

input_path = '/apps/gsc/matlab-library/scripts/depth_calibration/';
files = {'chaza-1_checkshot.txt','jodari-1_checkshot.txt','jodari-2_checkshot.txt','mdalasini-1_checkshot.txt','mzia-1_checkshot.txt','mzia-3_checkshot.txt'};
labels = strrep(files,'_',' ');
labels = strrep(labels,'.txt','');
seabed_d = [950 1150 1050 2280 1640 1780];
checkshot_xy = [667317.89,8873229.893;
                663648.27,8883553.664;
                660829.7,8881661.0;
                663887.50,8951659.50;
                660820.21,8907459.20;
                659941.26,8913297.58];

checkshot_datum= [0,-25,0,-25,0,0];

vel_meta_file = '/data/TZA/segy/2012_kussini_pgs_full_sequence/presdm/velocities/Kussini_Vint_depth_Isotropic/job_meta/';

%% read checkshot data and plot time-slowness for QC

num_checkshots = size(files,2);

max_t = 0;

figure;

for i = 1:num_checkshots
    checkshot_d_t{i} = dlmread([input_path files{i}]);
    checkshot_d_t{i}(:,1) = checkshot_d_t{i}(:,1)+checkshot_datum(i);
    checkshot_d_t{i}(:,2) = checkshot_d_t{i}(:,2)./2;         % convert from two-way time to one-way time
    checkshot_t_d{i} = fliplr(checkshot_d_t{i});              % create time-depth version to avoid confusion
    checkshot_t_sl{i} = dt2slowness(checkshot_d_t{i});
    
    scatter(checkshot_t_sl{i}(:,2),-checkshot_t_sl{i}(:,1));
    hold all;
    
    if checkshot_t_sl{i}(end,1) > max_t
        max_t = checkshot_t_sl{i}(end,1);
    end
    
%     scatter(checkshot_xy(i,1),checkshot_xy(i,2),'filled'); hold on;
    
end

% daspect([1 1 1]);

legend(labels,'Location','Best');

%% load velocity data

% read small cubes around each checkshot
% segy_read_binary to get trace definition
% use segy_index_byte_finder to get byte locations
% then read_traces_segy to get the traces

