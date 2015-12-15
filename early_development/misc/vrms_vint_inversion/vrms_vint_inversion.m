function vint = vrms_vint_inversion(aperture,n_samples,traces,tol,iter)
%final 
    % number of position to invert at once
    n_pos = (2*aperture+1)^2;

    % These parameters are the same for each group of traces
    % Create integration matrix
    C = cell(1,n_pos);
    [C{:}] = deal(sparse(tril(ones(n_samples,n_samples),0))^2);
    C = blkdiag(C{:});
    
    % Extract middle model trace from m
    if n_pos == 1
        start_crop = 1;
    else
        start_crop = 1+(median(1:1:n_pos)-1)*n_samples;
    end
        
    end_crop = start_crop + n_samples - 1;

    % Add lateral inversion contraint
    D = spdiags((-1)*ones(n_samples*n_pos,1),0,n_samples*n_pos,...
        n_samples*n_pos);
    E = repmat(spdiags(ones(n_samples,1),0,n_samples,...
        n_samples),n_pos,1);
    L = spalloc(n_samples*n_pos,n_samples*n_pos,...
        ((n_pos*2)-1)*n_samples)+D;
   
    L(:,start_crop:end_crop) = E; 
    
    G = C*L;
    
    %clear C D E L
    
    % can remove all matrices except G    

    %Process each group separately
    C_out = sparse(tril(ones(n_samples,n_samples)));
    
  
    
    ntraces = size(traces,2);

    % scale data appropriately
    t = (1:1:n_samples)';

    for ii=1:1:ntraces 
       traces(:,ii) = t.*(traces(:,ii).^2);
    end
    
    start_index = 1;
    n_groups = ntraces/n_pos;
    for ij=1:1:n_groups
        % fprintf('Group %d\n', ij);
        vrms = traces(:,start_index:n_pos*ij);
        vrms = vrms(:);
        [m,~] = lsqr(G,vrms,tol,iter);
        %clear vrms;

        m = C_out*m(start_crop:end_crop);
        vint(:,ij) = real(sqrt(m)); 

        %vint.pos(1,ij) = traces.pos(1,(start_index+(ij*n_pos))/2); 
        %vint.pos(2,ij) = traces.pos(2,(start_index+(ij*n_pos))/2);             
        start_index = start_index+n_pos;
    end

%     start_il = num2str(vint.pos(1,1));
%     end_il = num2str(vint.pos(1,end));
%     
%     start_xl = num2str(vint.pos(2,1));
%     end_xl = num2str(vint.pos(2,end));
%     
    % t_range
    
    % Write data
%     data_out = strcat(output_location,'/','vint_il_',start_il,'-',end_il,'_xl',start_xl,'-',end_xl);
%     file_vint = fopen(data_out,'a');
%     fwrite(file_vint,vint.m,'float32');
%     fclose(file_vint);

    % Write position information
    % meta_out = strcat(output_location,'/','vrms_pos_worker_',worker,'_block_',block);
    % dlmwrite(meta_out,vint.pos','delimiter','\t');

    %figure(i_block)
    %imagesc(vint.m)
    
end