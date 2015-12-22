
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

vel_files = {'chaza-1_depth_vint.txt','jodari-1_depth_vint.txt','jodari-2_depth_vint.txt','mdalasini-1_depth_vint.txt','mzia-1_depth_vint.txt','mzia-3_depth_vint.txt'};

datum = 5;

% figure;

for i = 1:num_checkshots
    model_d_v{i} = dlmread([input_path vel_files{i}]);
    model_d_v{i}(:,1) = model_d_v{i}(:,1)+datum;        % fix depths to start from 0
    model_d_v{i}(:,2) = model_d_v{i}(:,2)./1000;        % change units to m/ms
    model_t_d{i} = dv2tdepth(model_d_v{i});
    model_t_sl{i} = dv2slowness(model_d_v{i});
    
    seabed_t(i) = interp1(model_t_d{i}(:,2),model_t_d{i}(:,1),seabed_d(i));     % find time of seabed in model
    
%     scatter(model_t_d{i}(:,2),-model_t_d{i}(:,1));
%     hold all;
end

% Assume time-depth is correct down to seabed? 
% Or should water-column be first layer?
% How accurately are water depths measured at wells?

figure;

for i = 1:num_checkshots
    
    seabed_checkshot_t_d{i} = [checkshot_t_d{i}(:,1)-seabed_t(i),checkshot_t_d{i}(:,2)-seabed_d(i)];
    seabed_model_t_d{i} = [model_t_d{i}(:,1)-seabed_t(i),model_t_d{i}(:,2)-seabed_d(i)];
    
    seabed_model_t_error{i} = [seabed_checkshot_t_d{i}(:,1),interp1(seabed_model_t_d{i}(:,1),seabed_model_t_d{i}(:,2),seabed_checkshot_t_d{i}(:,1)) - seabed_checkshot_t_d{i}(:,2)];
    
    % scatter(seabed_model_t_sl{i}(:,2),-seabed_model_t_sl{i}(:,1),'filled'); hold on;
    scatter(seabed_model_t_error{i}(:,2),-seabed_model_t_error{i}(:,1),'filled'); hold on;
end

legend(labels,'Location','Best');

figure;

for i = 1:num_checkshots
    
    seabed_checkshot_t_sl{i} = [checkshot_t_sl{i}(:,1)-seabed_t(i),checkshot_t_sl{i}(:,2)];
    seabed_model_t_sl{i} = [model_t_sl{i}(:,1)-seabed_t(i),model_t_sl{i}(:,2)];
    
    seabed_model_t_sc{i} = [seabed_checkshot_t_sl{i}(:,1),seabed_checkshot_t_sl{i}(:,2)./(interp1(seabed_model_t_sl{i}(:,1),seabed_model_t_sl{i}(:,2),seabed_checkshot_t_sl{i}(:,1)))];
    
    seabed_model_sc_int{i} = interp1([seabed_model_t_sl{i}(1,1);seabed_model_t_sc{i}(:,1);seabed_model_t_sl{i}(end,1)],[seabed_model_t_sc{i}(1,2);seabed_model_t_sc{i}(:,2);seabed_model_t_sc{i}(end,2)],seabed_model_t_sl{i}(:,1));
    
    seabed_model_fix_t_sl = [seabed_model_t_sl{i}(:,1),seabed_model_t_sl{i}(:,2).*seabed_model_sc_int{i}];
    
    
    
    scatter(seabed_model_t_sc{i}(:,2),-seabed_model_t_sc{i}(:,1),'filled'); hold on;
    
end

legend(labels,'Location','Best');

% figure;
% 
% for i = 1:num_checkshots
%     
%        
%     abs_model_t_error{i} = [checkshot_t_d{i}(:,1),interp1(model_t_d{i}(:,1),model_t_d{i}(:,2),checkshot_t_d{i}(:,1)) - checkshot_t_d{i}(:,2)];
%     
%     scatter(abs_model_t_error{i}(:,2),-abs_model_t_error{i}(:,1)); hold all;
%     
% end


% figure;
% 
% for i = 1:num_checkshots
%     
%     checkshot_t_sl_int{i} = interp1(checkshot_t_sl{i}(:,1),checkshot_t_sl{i}(:,2),0:10:max_t,'linear','extrap');
%     model_t_sl_int{i} = interp1(model_t_sl{i}(:,1),model_t_sl{i}(:,2),0:10:max_t,'linear','extrap');
%     
%     
%     
%     error_t_sl{i} = model_t_sl_int{i} ./ checkshot_t_sl_int{i};
%     
%     plot(error_t_sl{i}); hold all;
% end
