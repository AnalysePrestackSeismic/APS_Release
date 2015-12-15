function sim_pre_stack_inversion(background_model,angle_stacks,...
    wavelets,angles,k,m,kc,mc,gamma,itermax,mu,tol,no_rho)

angle_stacks_path
n_angles
angles
wavelets

tol
iter
n_rho
c1
c2
c3
traces_process


% take average from well log if no background model

(filepaths,n_angles,aperture,n_samples,...
    n_blocks,byte_locations,block_proc_index,tol,iter,output_location)

% each worker processes n_block
 for i_block = 1:1:n_blocks
        %fprintf('\nBlock %d', i_block);
        % read in all the traces in the group
        start_index = block_proc_index(i_block,4)-block_proc_index(1,4)+1;
        end_index = start_index+block_proc_index(i_block,5)-block_proc_index(i_block,4);
        
        % this will need to loop over number of angle stacks
        angle_traces = segy_read_traces(filepath,n_samples,byte_locations,start_index,end_index);
        %traces = segy_read_traces(vrms_in,block_proc_index(i_block,4),block_proc_index(i_block,5));     
        ntraces = size(traces.data,2);

        
                               Lpmod = log(zpmod(:,xl));\n";
                      Lsmod = log(zsmod(:,xl))-(k*Lpmod)-kc;\n";
           
                                dLsmod(i,1) == -inf;
                                        dLsmod(i,1) = 0;
                                
             Ldmod = log(rhomod(:,xl))-(m*Lpmod)-mc;
                background_model = [Lpmod; dLsmod; dLdmod];
                  angle_stack_data(:,xl) = (1-mu)*angle_stack_data(:,xl) + mu*G*background_model;

        
        
        %fprintf('\nTotal traces %d \n', ntraces/n_pos);
[solution,flag,relres,itexit,resvec] = lsqr(G,angle_stack_data(:,xl),tol,itermax,[],[],background_model(:,xl));\n";

        start_index = 1;
        n_groups = ntraces/n_pos;
        for ij=1:1:n_groups
            % fprintf('Group %d\n', ij);
            vrms = traces.data(:,start_index:n_pos*ij);
            vrms = vrms(:);
            [m,~] = lsqr(G,vrms,tol,iter);
            %clear vrms;

            m = C_out*m(start_crop:end_crop);
            vint.m(:,ij) = real(sqrt(m)); 
            
            vint.pos(1,ij) = traces.pos(1,(start_index+(ij*n_pos))/2); 
            vint.pos(2,ij) = traces.pos(2,(start_index+(ij*n_pos))/2);             
            start_index = start_index+n_pos;
        end
        
        worker = num2str(block_proc_index(i_block,2));
        block = num2str(block_proc_index(i_block,3));
        
        % Write data
        data_out = strcat(output_location,'/','vrms_worker_',worker,'_block_',block);
        file_vint = fopen(data_out,'a');
        fwrite(file_vint,vint.m,'float32');
        fclose(file_vint);
        
        % Write position information
        meta_out = strcat(output_location,'/','vrms_pos_worker_',worker,'_block_',block);
        dlmwrite(meta_out,vint.pos','delimiter','\t');
        
        %figure(i_block)
        %imagesc(vint.m)
    
    end




% precondition data
angle_stacks = (1-mu)*angle_stacks + mu*G*background_model;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

