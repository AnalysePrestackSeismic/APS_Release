 for ii = 1:1:NumObjects
        if STATS(ii).Area > area_thresh_gt && STATS(ii).Area < area_thresh_lt
            fprintf('Creating surfaces for body %d of %d\n',ii,NumObjects);
            z_col = 2;
            skey_col = 1;
            pkey_col = 3;
            
            pos = STATS(ii).PixelList;
            sort_z_as = sortrows(pos,z_col);
            sort_z_ds = sortrows(pos,-z_col);
            pos_uniq = unique(pos(:,[pkey_col skey_col]),'rows');
            % col 1 skey
            % col 2 z
            % col 3 pekyn
            pos_body_surface_top(:,1) = pkeys(pos_uniq(:,1));
            pos_body_surface_top(:,2) = skeys(pos_uniq(:,2));
            
            [~,Locb] = ismember(pos_uniq,sort_z_ds(:,[pkey_col skey_col]),'rows');
            pos_body_surface_top(:,3) = (sort_z_ds(Locb,2)-1-10).*job_meta.s_rate/1000;
            
            pos_body_surface_bot(:,1) = pkeys(pos_uniq(:,1));
            pos_body_surface_bot(:,2) = skeys(pos_uniq(:,2));
            [~,Locb] = ismember(pos_uniq,sort_z_as(:,[3 1]),'rows');
            pos_body_surface_bot(:,3) = (sort_z_as(Locb,2)-1-10).*job_meta.s_rate/1000;
            z_mode = mode(pos_body_surface_top(:,3));
            % add mean depth to filename
            top_asc = strcat(out_dir_top,'/anom_surface_body_z',num2str(z_mode),'_area',num2str(STATS(ii).Area),'_id',num2str(ii),'_top.txt');
            bot_asc = strcat(out_dir_bot,'/anom_surface_body_z',num2str(z_mode),'_area',num2str(STATS(ii).Area),'_id',num2str(ii),'_bot.txt');
            
            area = repmat(STATS(ii).Area,size(pos_body_surface_bot(:,3),1),1);
            id = repmat(ii,size(pos_body_surface_bot(:,3),1),1);
            
            dlmwrite(top_asc,[pos_body_surface_top area id],'delimiter','\t','precision', 8)
            dlmwrite(bot_asc,[pos_body_surface_bot area id],'delimiter','\t','precision', 8)
            
            min_il = min(pos_body_surface_top(:,1));
            max_il = max(pos_body_surface_top(:,1));
            min_xl = min(pos_body_surface_top(:,2));
            max_xl = max(pos_body_surface_top(:,2));
            
            
            iii = 1;
            for il_ii = min_il:1:max_il
                poly_min(iii,:) = [il_ii,min(pos_body_surface_top((pos_body_surface_top(:,1) == il_ii),2))];
                poly_max(iii,:) = [il_ii,max(pos_body_surface_top((pos_body_surface_top(:,1) == il_ii),2))];
                iii = iii + 1;
            end
            
            for xl_ii = min_xl:1:max_xl
                poly_min(iii,:) = [min(pos_body_surface_top((pos_body_surface_top(:,2) == xl_ii),1)),xl_ii];
                poly_max(iii,:) = [max(pos_body_surface_top((pos_body_surface_top(:,2) == xl_ii),1)),xl_ii];
                iii = iii + 1;
            end
            
            poly = [poly_min;poly_max];
            poly = unique(poly,'rows');
            %poly = sortrows(poly);
            
            %poly_area = 
            %mid_il = min_il+floor((max_il-min_il)/2);
            %mid_xl = min_xl+floor((max_xl-min_xl)/2);
                        
            %[theta,~] = cart2pol(poly(:,1)-mid_il,poly(:,2)-mid_xl);
            %poly(:,3) = theta;
            %poly = sortrows(poly,3);
            poly(:,3) = pos_body_surface_top(ismember(pos_body_surface_top(:,1:2),poly(:,1:2),'rows'),3);
            pol_asc = strcat(out_dir_pol,'/anom_poly_z',num2str(z_mode),'_area',num2str(STATS(ii).Area),'_id',num2str(ii),'.txt');
            dlmwrite(pol_asc,poly,'delimiter','\t','precision', 8)         
           
            
            clear poly* pos sort_z_as sort_z_ds pos_uniq Locb pos_body_surface_top pos_body_surface_bot
            
        end
    end