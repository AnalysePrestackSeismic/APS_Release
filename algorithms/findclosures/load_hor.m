clear all
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
il_xl = 0;

start_point = pwd;  

% Load Horizons
for i_file = 1:1:nfiles  

hor_file = dlmread([filter_files.path{i_file},filter_files.names{i_file}]); % import of ascii horizon

if il_xl == 1
    min_il = min(hor_file(:,1));
    max_il = max(hor_file(:,1));
    min_xl = min(hor_file(:,2));
    max_xl = max(hor_file(:,2));

    il_inc = mode(diff(unique(hor_file(:,1))));                        % Primary Key (inline) increment (mode )
    xl_inc = mode(diff(unique(hor_file(:,2))));                        % Secondary Key (inline) Increment (mode)

    pkeyn = 1+(max_il-min_il)/il_inc;                 % Calculate Number of inlines (primary key)
    skeyn = 1+(max_xl-min_xl)/xl_inc;

    n_iline = (hor_file(:,1)-min_il)/il_inc+1;
    n_xline = (hor_file(:,2)-min_xl)/xl_inc+1;
    lin_ind = ((n_iline-1).*skeyn)+n_xline;

    slice = zeros(skeyn,pkeyn);
    slice(lin_ind) = hor_file(:,3);
else
    min_x = min(hor_file(:,1));
    max_x = max(hor_file(:,1));
    min_y = min(hor_file(:,2));
    max_y = max(hor_file(:,2));
    
    x_inc = mode(diff(unique(hor_file(:,1))));
    y_inc = mode(diff(unique(hor_file(:,2))));
    
    xkeyn = 1+(max_x-min_x)/x_inc;
    ykeyn = 1+(max_y-min_y)/y_inc;
    
    n_x = (hor_file(:,1)-min_x)/x_inc+1;
    n_y = (hor_file(:,2)-min_y)/y_inc+1;
    lin_ind = ((n_x-1).*ykeyn)+n_y;    
    
    slice = zeros(ykeyn,xkeyn);
    slice(lin_ind) = hor_file(:,3);   
    
    if i_file == 1;
        xs = min_x:x_inc:max_x;
        ys = min_y:y_inc:max_y;
    end
    
end

contint = 20;

[closures,closuregrd] = findclosures(slice,contint);

% Create contour display
figure(1)
map3=pmkmp(256);
subplot(1,3,1); imagesc(flipud(slice))

colormap(map3);  
colorbar
hold all
subplot(1,3,2); contour(slice,contint)

subplot(1,3,3); imagesc(flipud(closuregrd))
axis equal
axis tight

save_path = [filter_files.path{i_file},'results_contin_',num2str(contint),'_',filter_files.names{i_file},'/'];
mkdir(save_path);
saveas(1,[save_path,'img_',filter_files.names{i_file}], 'png');

%dir_poly = ['poly_',num2str(i_file)];
%mkdir(dir_poly);
%cd(dir_poly);
    %filter_files.names{i_file}]);
    
ii = 1;
for i_closure = 1:1:size(closures,1)   
    if size(closures{i_closure},2) > 40
        polycl = [xs(floor(closures{i_closure}(1,2:end)))',ys(floor(closures{i_closure}(2,2:end)))',repmat(closures{i_closure}(1,1),closures{i_closure}(2,1),1)];

%         if ii == 1;
%             figure(2)           
%             map3=pmkmp(256);
%             scatter(hor_file(:,1),hor_file(:,2),10,hor_file(:,3));
%             colormap(map3);  
%             colorbar
%             axis equal
%             hold all  
%             ii = 2;
%         else
%             figure(2)
%             scatter(polycl(:,1),polycl(:,2));    
%         end
        dlmwrite([save_path,'poly',num2str(i_closure),'-',num2str(closures{i_closure}(1)),'.txt'],polycl,'\t')
    end    
end
%saveas(2,[save_path,'img_poly',filter_files.names{i_file}], 'png');
%close all

cd(start_point);
fprintf('Completed horizon %d\n',i_file);
end