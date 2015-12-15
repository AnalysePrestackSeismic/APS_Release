function [] = minimum_energy_chi(block_mat_all,block_mat,process_files_mat,output_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    process_files = load(process_files_mat);
    [traces, process_positions] = segy_read_multiple_files(block_mat_all,process_files_mat, block_mat);

    for ii = 1:size(process_positions.ilxl_grid,1)
        if ii == 1
            [~, process_positions_ilxl_index_tmp] = ismember(process_positions.ilxl_grid(ii,1:2),flipud(process_positions.ilxl_pos(:,1:2)),'rows');
        else
            [~, process_positions_ilxl_index_tmp] = ismember(process_positions.ilxl_grid(ii,1:2),flipud(process_positions.ilxl_pos(process_positions_ilxl_index(ii-1,1)+process_positions.ilxl_step+1:end,1:2)),'rows');
        end
        process_positions_ilxl_index(ii,1) = length(process_positions.ilxl_pos)-process_positions_ilxl_index_tmp+1; 
    end
    
    ntraces_aperture = (1+2*process_positions.ilxl_aperture)^2;
    
    first_iter = 1;
    last_iter = length(process_positions_ilxl_index);
    
    for kk = first_iter:last_iter
        for ii = 1:process_files.nfiles
            traces_tmp = traces{ii};
            [ns,~] = size(traces_tmp);
            subvol_traces_tmp(1+(ii-1)*ns:ii*ns,:) = traces_tmp(...
                :,...
                process_positions_ilxl_index(kk)-((2*(process_positions.ilxl_aperture^2))+2*process_positions.ilxl_aperture):...
                process_positions_ilxl_index(kk)+((2*(process_positions.ilxl_aperture^2))+2*process_positions.ilxl_aperture));
        end
        data_tmp = reshape(subvol_traces_tmp,ns,[]);
        for jj = 1:ntraces_aperture
            data(:,1+(jj-1)*ns:jj*ns) = data_tmp(:,1+(jj-1)*process_files.nfiles:jj*process_files.nfiles)';
        end
        data(isnan(data)) = 0;
        NCava = [ones(process_files.nfiles,1) (sin(process_files.angle*pi/180).*sin(process_files.angle*pi/180))']\data;
        chi_raw = reshape((atand(-NCava(1,:)./NCava(2,:)))',ns,[]);
        live_traces_chi_raw = sum(isnan(chi_raw),1)~=ns;
        chi_raw = chi_raw(:,live_traces_chi_raw);
        chi_raw = [repmat((1:1:ns)',size(chi_raw,2),1),chi_raw(:)];
        chi_raw = chi_raw(~isnan(chi_raw(:,2)),:);
        chi_raw = chi_raw(and(chi_raw(:,2)<45,chi_raw(:,2)>0),:);
        chi_model(kk,1:2) = process_positions.ilxl_grid(kk,1:2);
        chi_model(kk,3:4) = ([ones(length(chi_raw),1) chi_raw(:,1)]\chi_raw(:,2))';
    end
    save(strcat(output_dir,sprintf('chi_model_block_%d',process_positions.block_id)),'chi_model');
end