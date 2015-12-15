start_dir = pwd;
load('int_grad_inv_proj_process_files.mat');
fid_ilxl = fopen('.tmp_ilxl','w');
fwrite(fid_ilxl,processing_grid.ilxl_grid(:,1:2),'int32');
fclose(fid_ilxl);
fid_ilxl = fopen('.tmp_ilxl','r');
single_ilxl = fread(fid_ilxl,'*single');
fclose(fid_ilxl);
single_ilxl = reshape(single_ilxl,[],2)';
total_cols = 0;
tmp_results = cell(4,1);

load(sprintf('int_grad_inv_proj_results_block_%d.mat',1));
cd('/data/KEN/dtect/kenya_3d_2012_L10ab_final/Seismics/binaries');
fid1 = fopen(strcat(results_out{1,1},'.bin'),'w');
fid2 = fopen(strcat(results_out{2,1},'.bin'),'w');
fid3 = fopen(strcat(results_out{3,1},'.bin'),'w');
fid4 = fopen(strcat(results_out{4,1},'.bin'),'w');
cd(start_dir);

for jj = 1:10
    total_cols_block = 0;
    fprintf('Reading blocks %d to %d into memory\n',1+(jj-1)*600,jj*600)
    if jj == 10
        for kk = 1:4
            tmp_results{kk,1} = single(zeros(2+n_samples_wb,599*2724+1104));
        end
    else
        for kk = 1:4
            tmp_results{kk,1} = single(zeros(2+n_samples_wb,600*2724));
        end
    end
    for ii = 1:600
        load(sprintf('int_grad_inv_proj_results_block_%d.mat',ii));
        [~,cols] = size(results_out{1,2});
        for kk = 1:4
            tmp_results{kk,1}(:,1+total_cols_block:total_cols_block+cols) = [single_ilxl(:,1+total_cols:total_cols+cols);single(results_out{kk,2})];
        end
        total_cols = total_cols + cols;
        total_cols_block = total_cols_block + cols;
        if ii/20 == floor(ii/20)
            fprintf('%d%% ',round(100*ii/600));
        end
    end
    cd('/data/KEN/dtect/kenya_3d_2012_L10ab_final/Seismics/binaries');
    fprintf('\nWriting blocks %d to %d to file\n',1+(jj-1)*600,jj*600)
    fwrite(fid1,tmp_results{1,1},'single');
    fwrite(fid2,tmp_results{2,1},'single');
    fwrite(fid3,tmp_results{3,1},'single');
    fwrite(fid4,tmp_results{4,1},'single');
    cd(start_dir);
end

fclose all;
