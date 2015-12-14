function [vol_3d] = load_3d_volume(job_meta_path,vol_index,start_slab,end_slab,aperture)
% function load 3D volume

job_meta = load(job_meta_path);                                 % Load job meta file

segy_plot_blocks(job_meta_path,'1');

pkey_inc_mode = mode(job_meta.pkey_inc);                        % Primary Key (inline) increment (mode )
skey_inc_mode = mode(job_meta.skey_inc);                        % Secondary Key (inline) Increment (mode)

pkeyn = 1+((job_meta.pkey_max(str2double(vol_index))-job_meta.pkey_min(str2double(vol_index)))...
    /job_meta.pkey_inc(str2double(vol_index)));                 % Calculate Number of inlines (primary key)
skeyn = 1+((job_meta.skey_max(str2double(vol_index))-job_meta.skey_min(str2double(vol_index)))...
    /job_meta.skey_inc(str2double(vol_index)));                 % Calculate Number of inlines (secondary key)

if str2double(end_slab) > job_meta.n_samples{str2double(vol_index)}
    end_slab = num2str(job_meta.n_samples{str2double(vol_index)});                    % update endslab if the data provided is of shorter length
end

n_slices = str2double(end_slab)-str2double(start_slab)+1;   % Number of Slices
vol_3d = zeros(pkeyn*skeyn,n_slices,'single');                                      % Initalize matrix for all inlines, xlines and slices

loopfin = size(job_meta.liveblocks,1);                          % Number of live blocks
lpi = 1;

while lpi <= loopfin
    i_block = job_meta.liveblocks(lpi);                            % Block Number for Current Live Block
    [~, traces, ilxl_read, ~] = ...
        node_segy_read(job_meta_path,vol_index,num2str(i_block));
    traces = [zeros(1,size(traces,2)); traces(2:end,:)];
     
    n_iline = (ilxl_read(:,1)-job_meta.pkey_min(str2double(vol_index)))/pkey_inc_mode+1;
    n_xline = (ilxl_read(:,2)-job_meta.skey_min(str2double(vol_index)))/skey_inc_mode+1;
    lin_ind = ((n_iline-1).*skeyn)+n_xline;    

    vol_3d(double(lin_ind),1:n_slices) = traces(str2double(start_slab):str2double(end_slab),:)';        
    
    fprintf('-- Block %d of %d --\n',lpi,loopfin);
    lpi = lpi + 1; 
end

vol_3d = reshape(vol_3d,skeyn,pkeyn,n_slices);
vol_3d = permute(vol_3d,[3 1 2]);

end

