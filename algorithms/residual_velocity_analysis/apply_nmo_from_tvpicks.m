function [] = apply_nmo_from_tvpicks(job_meta_path,velfile,i_block)
% -------------------------------------------------------------------------
%
% vel_file should be .mat file saved from residual velocity analysis
%
% this contains tvpairs cell matrix with picks


%
% to reduce printout in compilied version turned all warning off
warning off all;

%
% Load job meta information 
job_meta = load(job_meta_path);
%

% read all the data for this block
% node_segy_read(job_meta_path,vol_index,i_block)
[seismic, vol_traces, ilxl_read, offset_read] = node_segy_read(job_meta_path,'1',i_block);
% find the total number of offsets
offset = unique(offset_read);




% % reshape the gather array to the same 3d matrix as the angle volumes and
% drop as required
%fold = length(offset);
fold = max(seismic.trace_ilxl_bytes(:,7));
vol_traces = reshape(vol_traces,seismic.n_samples,fold,[]);
vol_out_traces = zeros(size(vol_traces),'single');


% figure(1); imagesc(squeeze(vol_traces(:,:,50))); colormap(gray); caxis([-2 2])

offset_read = reshape(offset_read,fold,{});


load(velfile);

% tnmo = 0.001*[vel_fn(:,1)];
% vnmo = [vel_fn(:,2)];
num_gathers = size(vol_traces,3);
stack = zeros(seismic.n_samples,num_gathers,'single');

%for gather = 1:num_gathers
for gather = 1:100
    
    input_gather = double(squeeze(vol_traces(:,:,gather)));
    
    inmo_gather = inmo(input_gather,0.002,double(offset_read(:,gather)),0,1500,1000);
    
    vol_out_traces(:,:,gather) = nmo(inmo_gather,0.002,double(offset_read(:,gather)),tvpairs2{gather}(:,1)./1000,tvpairs2{gather}(:,2).*1000,1000);
    
    stack(:,gather) = sum(vol_out_traces(:,:,gather),2);
    stack_fold = sum(vol_out_traces(:,:,gather) ~= 0,2);
    stack(:,gather) = stack(:,gather) ./ stack_fold;
    
    % figure(1); imagesc(input_gather); colormap(gray); caxis([-2 2])
    % figure(2); imagesc(vol_out_traces(:,:,gather)); colormap(gray); caxis([-2 2])
    
  
    
end

% for gather = 40:40:num_gathers;
%    plot_index = gather/40;
%    plot_traces(:,:,plot_index) = vol_out_traces(:,:,gather);
% end
for gather = 5:5:100;
   plot_index = gather/5;
   plot_traces(:,:,plot_index) = vol_out_traces(:,:,gather);
end

stack(isnan(stack)) = 0;
plot_traces(isnan(plot_traces)) = 0;
plot_traces = reshape(plot_traces,seismic.n_samples,fold*size(plot_traces,3));
figure; imagesc(plot_traces); colormap(gray); caxis([-0.5 0.5])
%figure; imagesc(stack); colormap(gray); caxis([-2 2])

vol_out_traces=vol_out_traces.*stack_mute;

clear plot_traces;


for gather = 5:5:100;
   plot_index = gather/5;
   plot_traces(:,:,plot_index) = vol_out_traces(:,:,gather);
end

stack(isnan(stack)) = 0;
plot_traces(isnan(plot_traces)) = 0;
plot_traces = reshape(plot_traces,seismic.n_samples,fold*size(plot_traces,3));
figure; imagesc(plot_traces); colormap(gray); caxis([-0.5 0.5])


% grab the actual angles from the gather to pick the correct indicies
% startvol = find(offset == startvol,1);
% endvol = find(offset == endvol,1);
% 
% % resize the ilxl data
% tmp_ilxlrd = ilxl_read(1:length(offset):end,:);
% ilxl_read = tmp_ilxlrd;
% clear tmp_ilxlrd;
% 
% %resize the offset header - no need as stacking
% %offset_read = offset_read(:,(startvol:volinc:endvol),:);
% 
% % now loop rjob_meta_pathound making however many angle gathers are requested
% aidx = 1;
% for kk = startvol:angwidth:endvol
%     % resize the traces data
%     vol_tracestmp = vol_traces(:,(startvol:volinc:(kk+angwidth)),:);
%     
%     for ii = 1:size(vol_tracestmp,3)
%         kds = vol_tracestmp(:,:,ii);
%         kdsb(:,ii) = sum((kds ~= 0),2);
%     end
%     % sum the gather and divide by live samples
%     %make logical of what is not zero and cumlatively sum it to get the fold
%     
%     angle_stk{aidx} = squeeze(sum(vol_tracestmp,2)) ./ kdsb;
%     aidx = aidx +1;
% end




% %% Save results
% aidx = 1;
% for kk = startvol:angwidth:endvol
%     
%     
%     resultno = 1;
%     % Save outputs into correct structure to be written to SEGY.
%     results_out{resultno,1} = 'Meta data for output files';
%     results_out{resultno,2}{1,1} = ilxl_read;
%     results_out{resultno,3} = 'is_gather'; % 1 is yes, 0 is no
%     %results_out{resultno,2}{2,1} = uint32(zeros(size(traces{vol_index_wb},2),1));
%     %was written as uint32(zeros(ntraces,1));
%     %results_out{resultno,2}{2,1} = offset_read';
%     
%     ebcstrtowrite = sprintf('%-3200.3200s',[results_out{resultno,1} '  ' ebdichdr '  ' tmpebc]);
%     results_out{resultno,1} = ebcstrtowrite;
%     
%     resultno = resultno + 1;
%     
%     
%     results_out{resultno,1} = strcat(testdiscpt);
%     %results_out{2,2} = digi_intercept;
%     results_out{resultno,2} = angle_stk{aidx};
%     results_out{resultno,3} = 0;
%     aidx = aidx +1;
%     
%     % check segy write functions - many different versions now!
%     if exist(strcat(job_meta.output_dir,'segy_io/'),'dir') == 0
%         output_dir = strcat(job_meta.output_dir,'segy_io/');
%         mkdir(output_dir);
%     else
%         output_dir = strcat(job_meta.output_dir,'segy_io/');
%     end
%     
%     i_block = str2double(i_block);
%     node_segy_write(results_out,i_block,job_meta.s_rate/1000,output_dir)
%     
% end

end


