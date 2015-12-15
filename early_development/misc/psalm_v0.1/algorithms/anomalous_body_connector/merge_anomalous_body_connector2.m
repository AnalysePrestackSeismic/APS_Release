function merge_anomalous_body_connector2(func_name,path_for_blocks)

start_point = pwd;
cd(path_for_blocks);  

process_files_mat = strcat(path_for_blocks,func_name,'_process_files.mat');
process_files = load(process_files_mat,'-mat');

fprintf('\n-- Joining results from %d blocks --\n',process_files.n_blocks); 
           
n_mat = 0;
while n_mat <= process_files.n_blocks
    system_for = sprintf('ls -B %stemp_%s_block* | wc -l',path_for_blocks,func_name);
    [~,n_mat] = system(system_for);  
    n_mat = str2double(n_mat);

    if n_mat == process_files.n_blocks % we have the files we need
        for i_block = 1:1:n_mat-1
            fprintf('-- Assessing connectivity at boundary between block %d and block %d --\n',... 
                i_block,i_block+1); 
                    % load result for block
                    if i_block == 1                
                        block_connect_mat = strcat(process_files.path_for_blocks,'temp_',...
                            process_files.func_name,'_block_',num2str(i_block),'.mat');
                        block_load{1} = load(block_connect_mat); % top block
                    else % make sure the updated values are taken
                        block_load{1} = block_load{2};
                        block_load{1}.I_id = reshape(block_load{1}.I_id,nz2,ntraces);
                        block_load{1}.I_vol = reshape(block_load{1}.I_vol,nz2,ntraces);
                        block_load{1}.I_depth = reshape(block_load{1}.I_depth,nz2,ntraces);
                        block_load{1}.I_azimuth_axis = reshape(block_load{1}.I_azimuth_axis,nz2,ntraces);
                        block_load{1}.I_azimuth_vals = reshape(block_load{1}.I_azimuth_vals,nz2,ntraces);
                    end

                    block_connect_mat = strcat(process_files.path_for_blocks,'temp_',...
                            process_files.func_name,'_block_',num2str(i_block+1),'.mat'); 
                    block_load{2} = load(block_connect_mat); % bottom block
                    
                    nz1 = size(block_load{1}.I_id,1);
                    nz2 = size(block_load{2}.I_id,1);
                    ntraces = size(block_load{1}.I_id,2)*size(block_load{1}.I_id,3);
                                     
                    % create boundary id matrix
                    % Column 1 = bottom boundary from top block
                    % Column 2 = top boundary from bottom block
                    
                    % Top block
                    connect_out.boundary_id(:,2*i_block-1) = reshape(block_load{1}.I_id(end,:),[],1);
                    connect_out.boundary_vol(:,2*i_block-1) = reshape(block_load{1}.I_vol(end,:),[],1);
                    connect_out.boundary_depth(:,2*i_block-1) = reshape(block_load{1}.I_depth(end,:),[],1);
                    connect_out.boundary_azimuth_vals(:,2*i_block-1) = reshape(block_load{1}.I_azimuth_vals(end,:),[],1);                    
                                        
                    % Bottom block
                    connect_out.boundary_id(:,2*i_block) = reshape(block_load{2}.I_id(1,:),[],1);
                    connect_out.boundary_vol(:,2*i_block) = reshape(block_load{2}.I_vol(1,:),[],1);
                    connect_out.boundary_depth(:,2*i_block) = reshape(block_load{2}.I_depth(1,:),[],1);
                    connect_out.boundary_azimuth_vals(:,2*i_block) = reshape(block_load{2}.I_azimuth_vals(1,:),[],1);
                    
                    % Find join locations across the boundary
                    logical_test = (connect_out.boundary_id(:,2*i_block-1).*connect_out.boundary_id(:,2*i_block) > 0);
                    ids(logical_test,1:2) = [connect_out.boundary_id(logical_test,2*i_block-1),...
                                            connect_out.boundary_id(logical_test,2*i_block)...
                                            ];
                    ids = unique(ids,'rows');
                    ids = ids(2:end,:);
                    
                    block_load{1}.I_id = block_load{1}.I_id(:);
                    block_load{2}.I_id = block_load{2}.I_id(:); 
                    
                    % Update ids in bottom block
                    % this will run out of memory so you need to block it
                    for i_id = 1:1:length(ids(:,1))
                        % Find index of ids that connect across the
                        % boundary in bottom block
                        positions_to_find = (block_load{2}.I_id == ids(i_id,2));
                        % Set the ids of these bodies to be the same as in
                        % top block
                        block_load{2}.I_id(positions_to_find) = ids(i_id,1);
                        
                        if (floor(i_id/5) == i_id/5)
                            fprintf('%d%% complete \n',round((i_id/length(ids(:,1))*100)));
                        elseif i_id == length(ids(:,1))
                            fprintf('Completed block %d. Saving Results... \n',i_block);
                        end
                    end
                                     
%                     block_sz = 50;
%                     h = zeros(length(block_load{1}.I_id),1);
%                     
%                     for ii=1:block_sz:length(ids(:,1))                        
%                         maxblock = ii + block_sz - 1;
%         
%                         if maxblock >  length(ids(:,1))
%                             maxblock = length(ids(:,1));
%                             block_sz = (maxblock - ii+1);    
%                         end
% 
%                         j = bsxfun(@eq,block_load{2}.I_id,ids(ii:maxblock,2)');    
%                         
%                         %f = j*diag(ids(ii:block_sz,1));
%                         f = bsxfun(@times,j,ids(ii:maxblock,1)'); % ids in bottom block now updated to ids in top block
%                         f = sum(f,2);
%                         g = block_load{2}.I_id.*bsxfun(@eq,sum(j,2),0); % find ids that need no updating
%                         h = h+g+f; % combine it all together 
%                         
%                         if (floor(ii/10) == ii/10)
%                             fprintf('%d%% complete \n',round((ii/length(ids(:,1))*100)));
%                         elseif ii == length(ids(:,1))
%                             fprintf('Completed block %d. Saving Results... \n',i_block);
%                         end
%     
%                     end
%                     
%                     block_load{2}.I_id = h;
%                     % Update boundary information to be consistent with
%                     % updates to ids
                    connect_out.boundary_id(:,2*i_block) = reshape(block_load{2}.I_id(1,:),[],1);
                    
                    clearvars ids 
                    % clearvars ids f j g h
                    
                    fprintf('- Saving results for boundary between block %d and block %d\n\n',... 
                        i_block,i_block+1);  
                    
                    %save block_load{1}
                    temp_results_mat = strcat(process_files.path_for_blocks,...
                            'temp_updated_',process_files.func_name,'_block_',...
                            num2str(i_block),'.mat');
                    struct_out.results_out{1,1} = 'body_id';
                    struct_out.results_out{1,2} = block_load{1}.I_id;
                    struct_out.anomalous_threshold = block_load{1}.anomalous_threshold;
                    save(temp_results_mat,'-struct','struct_out','-v7.3');
                    
                    if i_block == n_mat-1 % save final block
                        temp_results_mat = strcat(process_files.path_for_blocks,...
                            'temp_updated_',process_files.func_name,'_block_',...
                            num2str(i_block+1),'.mat');                           
                        struct_out.results_out{1,1} = 'body_id';
                        struct_out.results_out{1,2} = block_load{2}.I_id;
                        struct_out.anomalous_threshold = block_load{2}.anomalous_threshold;
                        save(temp_results_mat,'-struct','struct_out','-v7.3');
                    end                    
                     
        end        
   
        n_mat = process_files.n_blocks + 1; % break while loop 

        % Save boundary information
        connect_mat = strcat(process_files.path_for_blocks,...
            'connect_results.mat');
        save(connect_mat,'-struct','connect_out','-v7.3'); 

    end
end
cd(start_point);
end