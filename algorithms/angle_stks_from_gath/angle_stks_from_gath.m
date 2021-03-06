function [] = angle_stks_from_gath(job_meta_path,i_block,startvol,volinc,endvol,angwidth,tottracerun)
%
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
%% Parameters
% would normaly convert all parameters to double, but keep i_block as string as being passed to
% other modules; it does happen at the bottom of this program for output
%i_block = str2double(i_block);
%
% angle trace data ranges to use, vol is an angle trace either as a
% seperate angle volume or as a angle trace in an angle gather
% number of the first angle trace/volume to read
startvol = str2double(startvol);
% angle trace/volume increment
volinc = str2double(volinc);
% number of the last angle trace/volume to read
%endvol = job_meta.nvols;
endvol = str2double(endvol);
angwidth = str2double(angwidth);

% number of traces to run, put to zero to make it run all traces in the
% block, this is the default, this is also used to pass an inline (pkey
% number to use in testing has to be ilnnnn format
useselectemode = 0;

origstart_vol = startvol;
orig_endvol = endvol;
ebdichdr = ['segy io'];

if isempty(regexp(tottracerun,'il','once')) == 0
    useselectemode = 1;
    requiredinline =  str2double(regexprep(tottracerun,'il',''));
    %requiredinline =  str2double(strrep(tottracerun,'il',''));
    tottracerun = 0;
else
    tottracerun = str2double(tottracerun);
end
% tottracerun = 500;




% to reduce printout in compilied version turned all warning off
warning off all;

% end of parameters
%#####################################################################
%
% total number of volumes to load
totalvol = length(startvol:volinc:endvol);
droptraces = 1;
%
% Load job meta information 
job_meta = load(job_meta_path);
%
% if tottracerun == 0
%     ebdichdr = ['segy io'];
%     if useselectemode == 1;
%         testdiscpt = ['gathers_segy_io_',num2str(startvol),'_',num2str(volinc),'_',num2str(endvol)];
%     else
%         testdiscpt = ['gathers_segy_io_',num2str(startvol),'_',num2str(volinc),'_',num2str(endvol)];
%     end    
% else
%     ebdichdr = ['segy io'];
%     if useselectemode == 1;
%         testdiscpt = ['gathers_segy_io_',num2str(startvol),'_',num2str(volinc),'_',num2str(endvol)];
%     else
%         testdiscpt = ['gathers_segy_io_',num2str(startvol),'_',num2str(volinc),'_',num2str(endvol)];
%     end   
% end

% add the history of jobs run and this one to the curent ebcdic
if isfield(job_meta,'comm_history')
    ebdichdr2 = job_meta.comm_history;
    tmpebc = ebdichdr2{size(ebdichdr2,1),2};
else
    ebdichdr2{1,2} = '';
    tmpebc = '';
end

for ebcii = (size(ebdichdr2,1)-1):-1:1
    tmpebcc = regexp(ebdichdr2{ebcii,2},'/','split');
    tmpebc = [tmpebc tmpebcc{1}  tmpebcc{end}]; 
end
tmpebc = sprintf('%-3200.3200s',tmpebc);
clear tmpebcc ebdichdr2;
%
% Make ouput directories and create meta information
%fprintf('reading data for total volumes = %d\n',totalvol)
%
% read data to pick a water bottom on, this does not need to happen and then
% be ignored, should happen later after reading all the volumes
%

% read all the data for this block
% node_segy_read(job_meta_path,vol_index,i_block)
[~, vol_traces, ilxl_read, offset_read] = node_segy_read(job_meta_path,'1',i_block);
% find the total number of offsets
offset = unique(offset_read);


if droptraces == 1
    % % reshape the gather array to the same 3d matrix as the angle volumes and
    % drop as required
    nsamps = size(vol_traces,1);
    fold = length(offset);
    vol_traces = reshape(vol_traces,nsamps,fold,[]);
    % grab the actual angles from the gather to pick the correct indicies
    startvol = find(offset == startvol,1);
    endvol = find(offset == endvol,1);
    
    % resize the ilxl data
    tmp_ilxlrd = ilxl_read(1:length(offset):end,:);
    ilxl_read = tmp_ilxlrd;
    clear tmp_ilxlrd;
    
    %resize the offset header - no need as stacking
    %offset_read = offset_read(:,(startvol:volinc:endvol),:);
    
    % now loop round making however many angle gathers are requested
    aidx = 1;
    for kk = startvol:angwidth:endvol
        % resize the traces data
        vol_tracestmp = vol_traces(:,(kk:volinc:(kk+(angwidth-volinc))),:);
        kdsb = zeros(size(vol_tracestmp,1),size(vol_tracestmp,3));
        for ii = 1:size(vol_tracestmp,3)
            kds = vol_tracestmp(:,:,ii);
            kdsb(:,ii) = sum((kds ~= 0),2);
        end
        % sum the gather and divide by live samples
        %make logical of what is not zero and cumlatively sum it to get the fold
        
        angle_stk{aidx} = squeeze(sum(vol_tracestmp,2)) ./ kdsb;
        angle_stk{aidx}(isnan(angle_stk{aidx})) = 0;
        %figure(3); imagesc(angle_stk{aidx});  colormap(gray);
        aidx = aidx +1;
        clear kdsb
    end
end

i_block = str2double(i_block);

%% Save results
aidx = 1;
for kk = origstart_vol:angwidth:orig_endvol
    
    
    resultno = 1;
    % Save outputs into correct structure to be written to SEGY.
    results_out{resultno,1} = 'Meta data for output files';
    results_out{resultno,2}{1,1} = ilxl_read;
    results_out{resultno,3} = 'is_gather'; % 1 is yes, 0 is no
    results_out{resultno,2}{2,1} = uint32(zeros(size(angle_stk{aidx},2),1));
    %was written as uint32(zeros(ntraces,1));
    %results_out{resultno,2}{2,1} = offset_read';
     
    ebcstrtowrite = sprintf('%-3200.3200s',[results_out{resultno,1} '  ' ebdichdr '  ' tmpebc]);
    results_out{resultno,1} = ebcstrtowrite;
    
    resultno = resultno + 1;
    
    % correct file names added by SRW - 02/07/14
    
%     if kk == startvol
%         testdiscpt = ['gathers_segy_io_',num2str(origstart_vol),'_',num2str(volinc),'_',num2str(origstart_vol+angwidth)];
%     else
        testdiscpt = ['angle_stk_range_',num2str(kk),'_',num2str((kk+(angwidth-volinc)))];
%     end
       
    results_out{resultno,1} = strcat(testdiscpt);
    %results_out{2,2} = digi_intercept;
    results_out{resultno,2} = angle_stk{aidx};
    results_out{resultno,3} = 0;
    aidx = aidx +1;
    
    % check segy write functions - many different versions now!
    if exist(strcat(job_meta.output_dir,'bg_angle_stks/'),'dir') == 0
        output_dir = strcat(job_meta.output_dir,'bg_angle_stks/');
        mkdir(output_dir);
    else
        output_dir = strcat(job_meta.output_dir,'bg_angle_stks/');
    end
    
    
    node_segy_write(results_out,i_block,job_meta.s_rate/1000,output_dir)
    
end

end


