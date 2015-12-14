[vels_meta ilxl inp_vels] = segy_to_mat('189','193','/data/URY/segy/2013_pgs_uruguay_processing/BG_Total_Merge_Working/BG_Total_Velocity_Models_Raw_Merge_Limit2.segy');

% test triangulation

% everything has to be doubles for griddata

finl = double((floor(min(ilxl(:,1)))/5)*5);
linl = double((floor(max(ilxl(:,1)))/5)*5);

fxl = double((ceil(min(ilxl(:,2)))/6)*6-2);
lxl = double((ceil(max(ilxl(:,2)))/6)*6-2);


ilxl_double = double(ilxl);

out_vels = zeros(600,size(finl:5:linl,2),size(fxl:6:lxl,2));



% slice=200;
% out_vels(slice,:,:) = griddata(ilxl_double(:,2),ilxl_double(:,1),double(inp_vels(slice,:)),fxl:6:lxl,(finl:5:linl)');

parfor slice=1:600
    slice
    out_vels(slice,:,:) = griddata(ilxl_double(:,2),ilxl_double(:,1),double(inp_vels(slice,:)),fxl:6:lxl,(finl:5:linl)');
end

nan_check = ~sum(isnan(out_vels),1);

% meshgrid to get inline/crossline values of valid traces

[xl,il] = meshgrid(fxl:6:lxl,finl:5:linl);

ilxl_out = [il(nan_check),xl(nan_check)];

% put the output vel traces in a 2 dimension matrix

traces = out_vels(:,nan_check);


% write the output using super-easy and intuitive node_segy_write



output_vol{1,1} = 'Meta data for output files';
output_vol{2,2} = traces; % store output traces in inline/xline order
output_vol{1,2}{1,1} = ilxl_out; % store inline/xline values
output_vol{1,2}{2,1} = uint32(zeros(size(ilxl_out(:,1)))); % what is this for? offsets?
output_vol{1,1} = strcat('merged velocity traces ',date);
output_vol{1,3} = 'is_gather'; % 1 is yes, 0 is no
output_vol{2,3}=0;
output_vol{2,1} = 'velocity_merge_zone_v2';

output_dir = '/data/URY/segy/2013_pgs_uruguay_processing/BG_Total_Merge_Working/';
% check to see if the directory exists
if exist(output_dir,'dir') == 0
    mkdir(output_dir);
end/data



node_segy_write(output_vol,'1',16,output_dir); % write the output