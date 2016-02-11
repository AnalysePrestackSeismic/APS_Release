function [] = segy_qc_ava(vol_traces,ava,ns,input_angles)
%% ------------------ Disclaimer  ------------------
% 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) makes no representation or warranty, express or implied, in 
% respect to the quality, accuracy or usefulness of this repository. The code
% is this repository is supplied with the explicit understanding and 
% agreement of recipient that any action taken or expenditure made by 
% recipient based on its examination, evaluation, interpretation or use is 
% at its own risk and responsibility.
% 
% No representation or warranty, express or implied, is or will be made in 
% relation to the accuracy or completeness of the information in this 
% repository and no responsibility or liability is or will be accepted by 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) in relation to it.
%% ------------------ License  ------------------ 
% GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
%% github
% https://github.com/AnalysePrestackSeismic/
%% ------------------ FUNCTION DEFINITION ---------------------------------

i_sample = 180;
%i_trace = 50;

int = ava(1:ns,:);
grad = ava(1+ns:end,:);

for i_trace = 1:10:100
    plot((sind(input_angles)).^2,vol_traces(:,i_sample,i_trace),...
        '--rs',...
        'LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',10)
    hold on

    r_est = int(i_sample,i_trace)+grad(i_sample,i_trace).*(sind(input_angles)).^2;
    
    plot((sind(input_angles)).^2,r_est,...
        '--bs',...
        'LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','g',...
        'MarkerSize',10)
    title(['Sample ',num2str(i_sample)])            
    xlabel('sin^2(\Theta)')
    ylabel('Amplitude')  
    hold off
end


% -------------------------------------------------------------------------
%%
% job_meta = load(job_meta_path);
% start_sample = floor(str2double(start_time_wb)/(job_meta.s_rate/1000));
% end_sample = floor(str2double(end_time_wb)/(job_meta.s_rate/1000));
% % Read traces for 2/3 max angle stack
% vol_index_wb = ceil(job_meta.nvols*0.6667);
% [~, traces{vol_index_wb}, ilxl_read{vol_index_wb}] = ...
%     node_segy_read(job_meta_path,num2str(vol_index_wb),i_block);
% if size(traces{vol_index_wb},1) == 1 && size(traces{vol_index_wb},2) == 1 
%     return
% end

% Pick water bottom
% [wb_idx] = water_bottom_picker(traces{vol_index_wb},0);
% wb_idx(wb_idx < 0) = 1;
% win_sub = bsxfun(@plus,wb_idx,(0:job_meta.n_samples{vol_index_wb}-max(wb_idx))');
% win_ind = bsxfun(@plus,win_sub,(0:job_meta.n_samples{vol_index_wb}:...
% job_meta.n_samples{vol_index_wb}*(size(traces{vol_index_wb},2)-1)));
% n_samples = end_sample-start_sample+1;
% n_traces = size(traces{vol_index_wb},2);
% gather = zeros(job_meta.nvols,n_samples,n_traces);
% % Load block for remaining angle stacks
% for i_vol = 1:1:job_meta.nvols 
%     if i_vol ~= str2double(vol_index_wb) % don't repeat load for previously read stack
%         % Read traces
%         [~, traces{i_vol}, ~, ~] = ...
%         node_segy_read(job_meta_path,num2str(i_vol),i_block);
%     end   
%     % Flatten traces to water bottom
%     traces{i_vol} = traces{i_vol}(win_ind);  
%     gather(i_vol,1:n_samples,1:n_traces) = traces{i_vol}(start_sample:end_sample,:);
%     %traces{i_vol} = traces{i_vol}(,:);
%     input_angles(i_vol) = (job_meta.angle{i_vol}(2)+job_meta.angle{i_vol}(1))/2;
% 
%     figure(1)
%     subplot(1,job_meta.nvols,i_vol); imagesc(traces{i_vol});   
%     hold all
%     plot(repmat(start_sample,1,n_traces))
%     plot(repmat(end_sample,1,n_traces))
%     hold off
% end

% Load DIGI result

%for i_sample = 1:1:n_samples
    figure(2)
    subplot(ceil(sqrt(n_samples)),ceil(sqrt(n_samples)),i_sample);
    for i_trace=1:str2double(trace_decimate):n_traces
        if i_trace == 1
            plot((sind(input_angles)).^2,vol_traces(:,i_sample,i_trace),...
                            '--rs',...
                            'LineWidth',2,...
                            'MarkerEdgeColor','k',...
                            'MarkerFaceColor','g',...
                            'MarkerSize',10)
        else
            hold on
            plot((sind(input_angles)).^2,gather(:,i_sample,i_trace),...
                '--rs',...
                'LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10)
            hold off
        end
    end
    title(['Sample ',num2str(i_sample)])            
    xlabel('sin^2(\Theta)')
    ylabel('Amplitude')  
    if i_sample == n_samples
        subplot(ceil(sqrt(n_samples)),ceil(sqrt(n_samples)),i_sample+1);
        imagesc(gather(:,:,i_trace)');
    end
%end
% Test unflatten
% Unflatten data
% [ns,ntraces] = size(traces{1});
% traces_unflat = [traces{1};zeros(job_meta.n_samples{vol_index_wb}-ns,ntraces)];
% for kk = 1:length(wb_idx)
%     traces_unflat(:,kk) = circshift(traces_unflat(:,kk),wb_idx(kk));
% end
%% Plot AVA