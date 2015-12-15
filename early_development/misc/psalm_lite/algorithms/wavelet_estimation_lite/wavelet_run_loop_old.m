function [] = wavelet_run_loop(seismic_mat_path,i_block,n_blocks,output_dir)

% Code to set up wavelet_estimation_lite to run locally over X angle stacks divided into N blocks

%% Variables to be set:
seismic_mat_path = '/data/Global/dtect/SRW_test_dataset/SEGY/matlab/05-10_angle_stack_SRW.mat_lite';
output_dir = '/data/Global/dtect/SRW_test_dataset/SEGY/matlab/';


angle_stacks{1,1} = '/data/Global/dtect/SRW_test_dataset/SEGY/05-10_angle_stack_SRW.sgy';
angle_stacks{2,1} = '/data/Global/dtect/SRW_test_dataset/SEGY/10-15_angle_stack_SRW.sgy';
angle_stacks{3,1} = '/data/Global/dtect/SRW_test_dataset/SEGY/15-20_angle_stack_SRW.sgy';
angle_stacks{4,1} = '/data/Global/dtect/SRW_test_dataset/SEGY/20-25_angle_stack_SRW.sgy';
angle_stacks{5,1} = '/data/Global/dtect/SRW_test_dataset/SEGY/25-30_angle_stack_SRW.sgy';
angle_stacks{6,1} = '/data/Global/dtect/SRW_test_dataset/SEGY/30-35_angle_stack_SRW.sgy';
angle_stacks{7,1} = '/data/Global/dtect/SRW_test_dataset/SEGY/35-40_angle_stack_SRW.sgy';
angle_stacks{8,1} = '/data/Global/dtect/SRW_test_dataset/SEGY/40-45_angle_stack_SRW.sgy';

%i_block = '1';
%n_blocks = '5';
%%


 for a = 1:size(angle_stacks(:,1))
%     tic
%     
%     if a==1,
%     seismic_filepath{1,1} = '/data/Global/dtect/SRW_test_dataset/SEGY/05-10_angle_stack_SRW.sgy'  
%     else seismic_filepath{a,1} = sprintf('/data/Global/dtect/SRW_test_dataset/SEGY/%d-%d_angle_stack_SRW.sgy',a*5,a*5+5)
%     end
%     
%     
%     
%     
%     filename_index{a,1} = regexp(seismic_filepath,'(\d{2}-\d{2})','match');
%     
 seismic_filepath = angle_stacks(a,1)
 seismic.filepath = seismic_filepath;
    
    
    for i_block = 1:n_blocks
        wavelet_estimation_lite_sw(seismic_mat_path,i_block,n_blocks,output_dir,seismic_filepath)
        
        
    end
        
        
        
        if exist(strcat(output_dir,sprintf('wavelet_block_%d.mat',i_block)),'file') ~= 0
            load(strcat(output_dir,sprintf('wavelet_block_%d.mat',i_block)));
            if i_block==1
                tmp_w = w;
            else
                tmp_w(2:end,:) = tmp_w(2:end,:)+w(2:end,:);
                if sum(logical(w(1,:))) > sum(logical(tmp_w(1,:)))
                    tmp_w(1,:) = w(1,:);
                end
            end
        end
        
    
    
    %% Average wavelets
    avg_w = [tmp_w(1,:); bsxfun(@rdivide,tmp_w(3:end,:),tmp_w(2,:))];
    avg_w = avg_w(:,logical(1-logical(sum(isnan(avg_w)))));
    
    ns_win = 101;
    
    avg_w = avg_w(:,logical(sum(avg_w(2:end,:))));
    avg_w_time = [avg_w(1,:); circshift(ifft(avg_w(2:end,:),'symmetric'),floor(ns_win/2))];
    
    %   Save average wavelets
    save(strcat(char(output_dir),'average_wavelet_freq_',filename_index{1}{1}),'avg_w','-v7.3');
    save(strcat(char(output_dir),'average_wavelet_time_',filename_index{1}{1}),'avg_w_time','-v7.3');
    
    %   Compile final wavelets into cell array to be used in IG inversion
    all_wavelets_freq{1,a} = avg_w;
    all_wavelets_time{1,a} = avg_w_time;


%   Save final wavelet files
save(strcat(char(output_dir),'all_wavelets_freq'),'all_wavelets_freq','-v7.3');
save(strcat(char(output_dir),'all_wavelets_time'),'all_wavelets_time','-v7.3');

time = toc/60
%% Plot Wavelets
figure(1)

for a = 1: size(angle_stacks,1)
    
    if a == 1
        load('/data/Global/dtect/SRW_test_dataset/SEGY/matlab/test/average_wavelet_time_05-10.mat');
        subplot(8,1,1);
        imagesc(bsxfun(@rdivide,avg_w_time(2:end,:),max(avg_w_time(2:end,:))));
    else
        load(sprintf('/data/Global/dtect/SRW_test_dataset/SEGY/matlab/test/average_wavelet_time_%d-%d',a*5,a*5+5));
        subplot(8,1,a);
        imagesc(bsxfun(@rdivide,avg_w_time(2:end,:),max(avg_w_time(2:end,:))));
    end
    
    
end

end




%% Average wavelets across blocks

%        for ii = 1:str2double(n_blocks)
%             if exist(strcat(output_dir,sprintf('wavelet_block_%d.mat',ii)),'file') ~= 0
%                 load(strcat(output_dir,sprintf('wavelet_block_%d.mat',ii)));
%                 if ii==1
%                     tmp_w = w;
%                 else
%                     tmp_w(2:end,:) = tmp_w(2:end,:)+w(2:end,:);
%                     if sum(logical(w(1,:))) > sum(logical(tmp_w(1,:)))
%                         tmp_w(1,:) = w(1,:);
%                     end
%                 end
%             end
%         end
%
%         avg_w = [tmp_w(1,:); bsxfun(@rdivide,tmp_w(3:end,:),tmp_w(2,:))];
%         avg_w = avg_w(:,logical(1-logical(sum(isnan(avg_w)))));
%
%
%         avg_w = avg_w(:,logical(sum(avg_w(2:end,:))));
%         avg_w_time = [avg_w(1,:); circshift(ifft(avg_w(2:end,:),'symmetric'),floor(ns_win/2))];
%
%
%         save(strcat(char(output_dir),'average_wavelet_freq_',filename_index2{1}{1}),'avg_w','-v7.3');
%         save(strcat(char(output_dir),'average_wavelet_time_',filename_index2{1}{1}),'avg_w_time','-v7.3');