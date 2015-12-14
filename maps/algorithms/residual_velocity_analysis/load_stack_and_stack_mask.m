
        
trace_count = 0;
stack_all = zeros(2501,total_traces);
stack_mask_all = zeros(2501,total_traces);

for iblock = 2:23
    matfile = strcat('/segy/URY/2014_BG_water_column_imaging/matlab/rms_picks_test_block',int2str(iblock),'.mat');
    load(matfile,'stack','stack_mask','ntraces');
    trace_count = trace_count + ntraces;
    stack_all(:,trace_count-ntraces+1:trace_count) = stack;
    stack_mask_all(:,trace_count-ntraces+1:trace_count) = stack_mask;
end
    
% create vectors for scatter plot

[stky stkx] = find(stack_mask_all);
stky = stky.*2;
% decimated version

[stky10 stkx10] = find(stack_mask_all(:,10:10:end));
stky10 = stky10.*2;


