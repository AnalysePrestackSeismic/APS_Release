function report_anomalous_body_connector(path_for_blocks,vol_percent) 
% Create graphs and reports 
    start_point = pwd;
    cd(path_for_blocks);  
    
    func_name = 'anomalous_body_connector';

    process_files_mat = strcat(path_for_blocks,func_name,'_process_files.mat');
    process_files = load(process_files_mat,'-mat');
    
%     index_wb_mat = strcat(path_for_blocks,process_files.name(1),'_index_wb.mat');
%     load(index_wb_mat,'-mat');

    i_body_cum = 1; % this is so we can unwrap bodies from different blocks
    for i_block = 1:1:process_files.n_blocks      
        results_mat = strcat(process_files.path_for_blocks,...
            process_files.func_name,'_results_block_',...
            num2str(i_block),'.mat');         
        load(results_mat,'-mat'); 

        process_positions_mat = strcat(process_files.path_for_blocks,...
            process_files.func_name,'_positions_block_',...
            num2str(i_block),'.mat');  
        process_positions = load(process_positions_mat,'-mat'); 

        fprintf('- Extracting bodies of interest\n'); 

        % Create inline and crossline grids to enable their extraction
        nil = process_files.n_iline;
        nxl = process_files.n_xline;
        nz = length(process_positions.z_grid);
        
        il_grid = repmat(reshape(process_files.processing_grid.ilxl_grid(:,1),...
            1,nxl*nil),nz,1);
        xl_grid = repmat(reshape(process_files.processing_grid.ilxl_grid(:,2),...
            1,nxl*nil),nz,1);
        
        vols_out = unique(results_out{2,2});
        n_vols = length(vols_out);
        vol_cut_off = vols_out(end-floor(n_vols*vol_percent/100));
        
        fprintf('- Applying a volume cut-off of %d%% (%d)\n',vol_percent,vol_cut_off); 
        
        positions_vol_cut = (results_out{2,2} > vol_cut_off);
        ids_out = unique(nonzeros(results_out{1,2}(positions_vol_cut)));
        for i_body = 1:1:length(ids_out) 
            % Get positions based on body id of interest
            positions_to_report = (results_out{1,2} == ids_out(i_body));    
            % Create inline and crossline grids to enable their extraction

            % Create time grid different for each block
            % ask for flattened times
            z_grid = repmat((process_positions.z_grid-1)*...
                process_files.s_rate/1000,1,nxl*nil); %+index_wb;

            % Body position information
            pos_body{i_body_cum} = [il_grid(positions_to_report) ...
            xl_grid(positions_to_report) ...
            z_grid(positions_to_report) ...
            results_out{1,2}(positions_to_report) ... % id
            results_out{2,2}(positions_to_report) ... % volume
            results_out{3,2}(positions_to_report) ... % depth
            results_out{4,2}(positions_to_report)]; % azimuth

            i_body_cum = i_body_cum + 1;
        end 
    end
        
    fprintf('- Extracting surfaces from bodies of interest\n'); 
    % Combine extracted body information
    pos_body = cell2mat(pos_body'); 
    pos_body = sortrows(pos_body,4); % sort on id
    
    % save bodies
    pos_body_mat = strcat(process_files.path_for_blocks,process_files.func_name,...
        '_volume_cut-off_bodies_results.mat'); 
    save(pos_body_mat,'pos_body','-v7.3');    
    
    ids_out = unique(pos_body(:,4)); % get ids 
    report_out = zeros(length(ids_out),7);

    for i_body = 1:1:length(ids_out)
        fprintf('-- Creating statistics and contouring for body %d of %d --\n',i_body,length(ids_out)); 
        % create stats for report
        % extract single body
        body_pos = pos_body(pos_body(:,4) == ids_out(i_body),:);
        il_out = median(body_pos(:,1));
        xl_out = median(body_pos(:,2));
        report_out(i_body,1:6) = [body_pos(1,4)... % id
                il_out(1) ... % middle inline
                xl_out(1) ... % middle crossline
                body_pos(1,5) ... % volume
                body_pos(1,6) ... % depth
                body_pos(1,7)]; % azimuth

        % extract surfaces    
        sort_z_as = sortrows(body_pos,3);
        sort_z_ds = sortrows(body_pos,-3);
        pos_uniq = unique(body_pos(:,1:2),'rows');

        [~,Locb] = ismember(pos_uniq,sort_z_ds(:,1:2),'rows');         
        pos_body_surface{i_body,1} = [pos_uniq sort_z_ds(Locb,3)]; % top

        [~,Locb] = ismember(pos_uniq,sort_z_as(:,1:2),'rows');                          
        pos_body_surface{i_body,2} = [pos_uniq sort_z_as(Locb,3)]; % base

        % Save surfaces
        top_asc = strcat(process_files.path_for_blocks,'/surface_body_',num2str(i_body),'top.txt'); 
        bot_asc = strcat(process_files.path_for_blocks,'/surface_body_',num2str(i_body),'bot.txt'); 
        dlmwrite(top_asc,pos_body_surface{i_body,1},'delimiter','\t','precision', 8)
        dlmwrite(bot_asc,pos_body_surface{i_body,2},'delimiter','\t','precision', 8)

        % Contour tops
        x = unique(pos_body_surface{i_body,1}(:,1)); % inline                                
        y = unique(pos_body_surface{i_body,1}(:,2)); % crossline                              
        
        if length(x) > 2 && length(y) > 2
            [X,Y] = meshgrid(x,y);                            
            [Lia2,Locb2] = ismember([X(:) Y(:)],pos_body_surface{i_body,1}(:,1:2),'rows');                          
            z = NaN(length(y)*length(x),1);
            z(Lia2) = pos_body_surface{i_body,1}(nonzeros(Locb2),3);
            z = reshape(z,length(y),length(x));                            
            pos_body_contour{i_body,1} = contourc(z);


            h1a = imagesc(x,y,z);
            hold all
            contour(x,y,z,'--w','LineWidth',2) 
            hold off
            title(sprintf('Body: %d, Input volume: %s,\nAnomalous threshold: %f, \n Volume cut-off: %f' ...
                ,i_body,process_files.name{1},anomalous_threshold,vol_cut_off),'interpreter','none');
            xlabel('Inline')
            ylabel('Crossline')
            colormap(flipud(jet))
            c = colorbar;
            ylabel(c,'Depth') 
            %text(6500,6500,report_header)
            contour_fig = strcat('contour_body_',num2str(i_body),'.pdf'); 
            saveas(h1a,contour_fig,'pdf') 
      

            % Test for closing contour
            jump = 1;
            nclosures = 0;                            
            while jump < length(pos_body_contour{i_body,1})                                
                if pos_body_contour{i_body,1}(1,jump+1) == ...
                        pos_body_contour{i_body,1}(1,jump+pos_body_contour{i_body,1}(2,jump)) ...
                        && pos_body_contour{i_body,1}(2,jump+1) == ...
                        pos_body_contour{i_body,1}(2,jump+pos_body_contour{i_body,1}(2,jump))
                    nclosures = nclosures+1;                                    
                end
                report_out(i_body,7) = nclosures;
                % add lowest closure contour
                jump = jump+pos_body_contour{i_body,1}(2,jump)+1;              
            end
        end
    end

    % sort on volume
    fprintf('\n-- Creating report and graphs --\n'); 
    report_out = sortrows(report_out,-4);
    mode_depth = mode(report_out(:,5));    
    mode_azimuth = mode(report_out(:,6)); 
    max_vol = max(report_out(:,4));
    min_vol = min(report_out(:,4));

    name = char(process_files.name(1));
    report_header{1,:} = sprintf('Input volume: %s',name); 
    report_header{2,:} = sprintf('Anomalous threshold: %f', anomalous_threshold);
    report_header{3,:} = sprintf('Modal depth: %f', mode_depth);
    report_header{4,:} = sprintf('Modal azimuth: %f', mode_azimuth);
    report_header{5,:} = sprintf('Min volume: %f', min_vol);
    report_header{6,:} = sprintf('Max volume: %f', max_vol);
    report_header{7,:} = sprintf('Volume cut off: %f',vol_cut_off);
    report_header{8,:} = sprintf('Body ID\tInline\tCrossline\tVolume\tDepth\tAzimuth\tNumber of closures');

%     h1 = scatter(report_out(:,2),report_out(:,3),report_out(:,4),report_out(:,5),'filled'); 
%     title(sprintf('Input volume: %s,\nAnomalous threshold: %f, \n Volume cut-off: %f' ...
%     ,process_files.name{1},anomalous_threshold,vol_cut_off),'interpreter','none');
%     xlabel('Inline')
%     ylabel('Crossline')
%     colormap(flipud(jet))
%     c = colorbar;
%     ylabel(c,'Depth') 
%     %text(6500,6500,report_header)
%     saveas(h1,'scatter_blobs_bodies.pdf','pdf') 

    % pos_body{i_body_cum} = [il_grid(positions_to_report) ...
    %     xl_grid(positions_to_report) ...
    %     z_grid(positions_to_report) ...
    %     results_out{1}(positions_to_report) ... % id
    %     results_out{2}(positions_to_report) ... % volume
    %     results_out{3}(positions_to_report) ... % depth
    %     results_out{4}(positions_to_report)]; % azimuth                  

%     h2 = scatter(pos_body(:,1),pos_body(:,2),10,pos_body(:,6),'filled'); 
%     title(sprintf('Input volume: %s,\nAnomalous threshold: %f, \n Volume cut-off: %f' ...
%     ,process_files.name{1},anomalous_threshold,vol_cut_off),'interpreter','none');
%     xlabel('Inline')
%     ylabel('Crossline')
%     colormap(flipud(jet))
%     c = colorbar;
%     ylabel(c,'Depth')
%     saveas(h2,'scatter_bodies_mode_depth.pdf','pdf') 

%     h3 = scatter(pos_body(:,1),pos_body(:,2),10,pos_body(:,5),'filled'); 
%     title(sprintf('Input volume: %s,\nAnomalous threshold: %f, \n Volume cut-off: %f' ...
%     ,process_files.name{1},anomalous_threshold,vol_cut_off),'interpreter','none');
%     xlabel('Inline')
%     ylabel('Crossline')
%     colormap(flipud(jet))
%     c = colorbar;
%     ylabel(c,'Volume')
%     %text(6500,6500,report_header)
%     saveas(h3,'scatter_bodies_mode_volume.pdf','pdf') 

%     h4 = scatter3(pos_body(:,1),pos_body(:,2),pos_body(:,3),25,pos_body(:,5),'filled');
%     title(sprintf('Input volume: %s,\nAnomalous threshold: %f, \n Volume cut-off: %f' ...
%     ,process_files.name{1},anomalous_threshold,vol_cut_off),'interpreter','none');
%     xlabel('Inline')
%     ylabel('Crossline')
%     zlabel('Z')
%     colormap(flipud(jet))
%     c = colorbar;
%     ylabel(c,'Volume')
%     set(gca,'ZDir','rev')
%     saveas(h4,'scatter_bodies_3d.pdf','pdf')

    report_mat = strcat(process_files.path_for_blocks,'report_connected_bodies.mat'); 
    save(report_mat,'report_header','report_out','-v7.3');

    report_asc = strcat(process_files.path_for_blocks,'connected_bodies.txt'); 
    dlmwrite(report_asc,report_header{1,:},'delimiter','')
    for i_out = 2:1:length(report_header)
        dlmwrite(report_asc, report_header{i_out,:},'delimiter','','-append')
    end
    dlmwrite(report_asc, report_out,'delimiter','\t','precision', 8,'-append')
    % create graphs 
    
cd(start_point);
end
