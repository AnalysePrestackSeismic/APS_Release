function [] = apply_nmo_single_function(job_meta_path,vel_file,i_block)
% -------------------------------------------------------------------------
% Authors: Charles Jones 
% -------------------------------------------------------------------------
%% Parameters
% would normaly convert all parameters to double, but keep i_block as string as being passed to
% other modules; it does happen at the bottom of this program for output
%i_block = str2double(i_block);
%
% angle trace data ranges to use, vol is an angle trace either as a
% seperate angle volume or as a angle trace in an angle gather
% number of the first angle trace/volume to read
% startvol = str2double(startvol);
% % angle trace/volume increment
% volinc = str2double(volinc);
% % number of the last angle trace/volume to read
% %endvol = job_meta.nvols;
% endvol = str2double(endvol);
% angwidth = str2double(angwidth);

% % number of traces to run, put to zero to make it run all traces in the
% % block, this is the default, this is also used to pass an inline (pkey
% % number to use in testing has to be ilnnnn format
% useselectemode = 0;
% 
% if isempty(regexp(tottracerun,'il','once')) == 0
%     useselectemode = 1;
%     requiredinline =  str2double(regexprep(tottracerun,'il',''));
%     %requiredinline =  str2double(strrep(tottracerun,'il',''));
%     tottracerun = 0;
% else
%     tottracerun = str2double(tottracerun);
% end
% % tottracerun = 500;




% to reduce printout in compilied version turned all warning off
warning off all;

% end of parameters
%#####################################################################
%
% % total number of volumes to load
% totalvol = length(startvol:volinc:endvol);
% droptraces = 1;
%
% Load job meta information 
job_meta = load(job_meta_path);
%
% if tottracerun == 0
%     ebdichdr = ['segy io'];
%     if useselectemode == 1;
%         testdiscpt = ['gathers_segy_io_',num2str(startvol),'_',num2str(volinc),'_',num2str(endvol)];
%     else
%         testdiscpt = ['gathers_segy_io_',num2str(startvol),'_',num2str(volinc),'i_block_',num2str(endvol)];
%     end    
% else
%     ebdichdr = ['segy io'];
%     if useselectemode == 1;
%         testdiscpt = ['gathers_segy_io_',num2str(startvol),'_',num2str(volinc),'_',num2str(endvol)];
%     else
%         testdiscpt = ['gathers_segy_io_',num2str(startvol),'_',num2str(volinc),'_',num2str(endvol)];
%     end   
% end

% % add the history of jobs run and this one to the curent ebcdic
% if isfield(job_meta,'comm_history')
%     ebdichdr2 = job_meta.comm_history;
%     tmpebc = ebdichdr2{size(ebdichdr2,1),2};
% else
%     ebdichdr2{1,2} = '';
%     tmpebc = '';
% end
% 
% for ebcii = (size(ebdichdr2,1)-1):-1:1
%     tmpebcc = regexp(ebdichdr2{ebcii,2},'/','split');
%     tmpebc = [tmpebc tmpebcc{1}  tmpebcc{end}]; 
% end
% tmpebc = sprintottraceruntf('%-3200.3200s',tmpebc);
% clear tmpebcc ebdichdr2;
%
% % Make ouput directories and create meta information
% fprintf('reading data for total volumes = %d\n',totalvol)
%
% read data to pick a water bottom on, this does not need to happen and then
% be ignored, should happen later after reading all the volumes
%

% read all the data for this block
% node_segy_read(job_meta_path,vol_index,i_block)
[seismic, vol_traces, ilxl_read, offset_read] = node_segy_read(job_meta_path,'1',i_block);
% find the total number of offsets
offset = unique(offset_read);



% % reshape the gather array to the same 3d matrix as the angle volumes and
% drop as required
nsamps = size(vol_traces,1);
%fold = length(offset);
fold = max(seismic.trace_ilxl_bytes(:,7));
vol_traces = reshape(vol_traces,nsamps,fold,[]);
vol_out_traces = zeros(size(vol_traces),'single');


% figure(1); imagesc(squeeze(vol_traces(:,:,50))); colormap(gray); caxis([-2 2])

offset_read = reshape(offset_read,fold,{});

vel_fn = csvread(vel_file);


tnmo = 0.001*[vel_fn(:,1)];
vnmo = [vel_fn(:,2)];
num_gathers = size(vol_traces,3);
stack = zeros(nsamps,num_gathers,'single');

for gather = 1:num_gathers
    
    input_gather = double(squeeze(vol_traces(:,:,gather)));
    
    vol_out_traces(:,:,gather) = nmo(input_gather,0.002,double(offset_read(:,gather)),tnmo,vnmo,40);
    
    stack(:,gather) = sum(vol_out_traces(:,:,gather),2);
    stack_fold = sum(vol_out_traces(:,:,gather) ~= 0,2);
    stack(:,gather) = stack(:,gather) ./ stack_fold;
    
    % figure(1); imagesc(input_gather); colormap(gray); caxis([-2 2])
    % figure(2); imagesc(vol_out_traces(:,:,gather)); colormap(gray); caxis([-2 2])
    
  
    
end

for gather = 40:40:num_gathers;
   plot_index = gather/40;
   plot_traces(:,:,plot_index) = vol_out_traces(:,:,gather);
end

plot_traces = reshape(plot_traces,nsamps,fold*size(plot_traces,3));
figure(2); imagesc(plot_traces); colormap(gray); caxis([-2 2])
figure(3); imagesc(stack); colormap(gray); caxis([-2 2])

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


