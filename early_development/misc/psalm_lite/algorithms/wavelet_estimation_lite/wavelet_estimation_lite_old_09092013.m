
function [] = wavelet_estimation_lite(seismic_mat_path,i_block,n_blocks,output_dir)

%%
% INPUT ARGUEMENTS
% Seismic_mat_path = Geometry file containing In Xline Byte in Binary
% Format (as output by segy_make_structure_lite)
% i_block = ith block of total n blocks.
% n_block = n blocks of seismic volume
% output_dir = directory to write output to.
%% READING IN SEISMIC VOLUME

%   Seismic = load(seismic_mat_path);     

    seismic.filepath = '/data/TZA/dtect/Ed_Cliffe_MSc/matlab/inputs/kusini_2012_113j05_40-45_angle_stk_merge.segy';

    tmp_seismic = fread(fopen(seismic_mat_path,'r'),'double');
    
%   Collect number of samples, ILXL bytes and number of traces

    seismic.n_samples = tmp_seismic(1,1);
    seismic.trace_ilxl_bytes = reshape(tmp_seismic(2:end),3,[])';
    seismic.n_traces = size(seismic.trace_ilxl_bytes,1);

    if ischar(i_block)
        i_block = str2double(i_block);
    end
    if ischar(n_blocks)
        n_blocks = str2double(n_blocks);
    end

    dec = 4;
    
    %Read traces
    [traces, ilxl_read] = node_segy_read_traces_lite(seismic,i_block,n_blocks,dec);
    
    load('/data/TZA/dtect/Ed_Cliffe_MSc/matlab/inputs/wb.mat');
    wb_idx = double(wb_idx(ismember(wb_idx(:,1:2),ilxl_read,'rows'),3))';
%     wb_idx = water_bottom_flatten_lite(traces);

    n_traces_in_block = size(traces,2);

    ns_win = 101;
    ns_overlap = 50;
    n_win = 1+floor((seismic.n_samples-ns_win)/(ns_win-ns_overlap));
    w = zeros(2+ns_win,n_win);

    win_idx = (0:1:ns_win-1)';

    for ii = 1:n_win
        win_sub = bsxfun(@plus,wb_idx,win_idx+(ii-1)*ns_overlap);    
        win_ind = bsxfun(@plus,win_sub,(0:seismic.n_samples:seismic.n_samples*(n_traces_in_block-1)));
        win_ind = win_ind(:,sum(win_sub > seismic.n_samples)==0);
        if isempty(win_ind)
            break
        end
        w(1,ii) = floor(ns_win/2)+(ii-1)*ns_overlap;
        w(2,ii) = size(win_ind,2);
        w(3:end,ii) = sum(abs(fft(traces(win_ind))),2);
    end

    save(strcat(output_dir,'wavelet_block_',num2str(i_block)),'w','-v7.3');
    
    if i_block == n_blocks 
        system_call = strcat(sprintf('ls -B  %s',output_dir),'wavelet_block_*.mat | wc -l');
        [~,result_count] = system(system_call);
        result_count = str2double(result_count);
        
        while result_count ~= n_blocks
            system_call = strcat('ls -B ',output_dir,'wavelet_block_*.mat | wc -l');
            [~,result_count] = system(system_call);
            result_count = str2double(result_count);
        end
        
       for ii = 1:n_blocks
            if exist(strcat(output_dir,sprintf('wavelet_block_%d.mat',ii)),'file') ~= 0
                load(strcat(output_dir,sprintf('wavelet_block_%d.mat',ii)));
                if ii==1
                    tmp_w = w;
                else
                    tmp_w(2:end,:) = tmp_w(2:end,:)+w(2:end,:);
                    if sum(logical(w(1,:))) > sum(logical(tmp_w(1,:)))
                        tmp_w(1,:) = w(1,:);
                    end
                end
            end
        end
        
        avg_w = [tmp_w(1,:); bsxfun(@rdivide,tmp_w(3:end,:),tmp_w(2,:))];
        avg_w = avg_w(:,logical(1-logical(sum(isnan(avg_w)))));
        avg_w = avg_w(:,logical(sum(avg_w(2:end,:))));

        avg_w_time = [avg_w(1,:); circshift(ifft(avg_w(2:end,:),'symmetric'),floor(ns_win/2))];

        save(strcat(output_dir,'average_wavelet_freq'),'avg_w','-v7.3');  
        save(strcat(output_dir,'average_wavelet_time'),'avg_w_time','-v7.3');
    end
end