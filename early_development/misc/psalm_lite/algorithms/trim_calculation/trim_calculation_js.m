function [] = trim_calculation_js(seismic_mat_path,i_block,n_blocks,seis_ext,output_dir,sample_rate,shift,zsmooth)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% trim_calculation('/data/NOR/segy/hegg_fwi_test_2012/gathers/run3_hegg_fwi_model_vti_full_area_starting_model.mat_lite','/data/NOR/segy/hegg_fwi_test_2012/gathers/run3_hegg_fwi_model_vti_full_area_starting_model.segy',100,1000000,0,'',12,50)

    i_block = str2double(i_block);
    n_blocks = str2double(n_blocks);
    sample_rate = str2double(sample_rate);
    shift = str2double(shift);
    zsmooth = str2double(zsmooth);
    
    % Open the scan segy file
    %tmp_seismic = fread(fopen(seismic_mat_path,'r'),'double');
    fid = fopen(seismic_mat_path,'r');
    seismic.n_samples = fread(fid,[1,1],'double');
    seismic.trace_ilxl_bytes = fread(fid,[3,inf],'double'); %reshape(tmp_seismic(2:end),3,[])';
    seismic.trace_ilxl_bytes = seismic.trace_ilxl_bytes';
    seismic.n_traces = size(seismic.trace_ilxl_bytes,1);
    
    if isempty(seis_ext) == 1
        seismic.filepath = regexprep(seismic_mat_path, '.mat_lite', seis_ext);
    else
        seismic.filepath = regexprep(seismic_mat_path, 'mat_lite', seis_ext);     
        seis_ext = strcat('.',seis_ext);
    end
    
    % Read some data
    dec = 0;
    [traces, ilxl_read, offset_read] = node_segy_read_traces_lite_gathers(seismic,i_block,n_blocks,dec);
    [n_samples,n_traces] = size(traces);
    
    % Make results meta information
    results_out{1,1} = 'Meta data for output files';
    results_out{1,2}{1,1} =  ilxl_read; %[]
    results_out{1,2}{2,1} =  offset_read;
    results_out{1,2}{1,2} =  ilxl_read;
    results_out{1,2}{2,2} =  offset_read;
    
    filename = regexprep(seismic.filepath, seis_ext, '');
    indexfile = strfind(filename, '/');
    filename = filename(indexfile(end)+1:end);
      
    results_out{2,1} = strcat(filename,'_trim_shifts_',num2str(shift),'-',num2str(zsmooth));
    results_out{3,1} = strcat(filename,'_trim_data_',num2str(shift),'-',num2str(zsmooth));
    results_out{4,1} = strcat(filename,'_trim_sum_',num2str(shift),'-',num2str(zsmooth));
    results_out{2,2} = zeros(n_samples,n_traces);
    results_out{3,2} = zeros(n_samples,n_traces);
     
    % loop through all traces and select gathers
    start_idx = 1;
    end_idx = start_idx;
    trace_counter = 1;    
    while start_idx < n_traces
        if sum(ilxl_read(start_idx,:) == ilxl_read(end_idx,:)) == 2 && end_idx < n_traces;   
            end_idx = end_idx + 1;
        elseif end_idx-start_idx == 1
            % test for single fold gather
            results_out{2,2}(:,start_idx:end_idx-1) = zeros(n_samples,1);
            results_out{3,2}(:,start_idx:end_idx-1) = zeros(n_samples,1);
            results_out{4,2}(:,trace_counter) = zeros(n_samples,1);
            trace_counter = trace_counter + 1;
            start_idx = start_idx + 1;
            end_idx = start_idx;
        else
        end_idx = end_idx - 1;
        % Extract inline and crossline information    
        results_out{1,2}{1,3}(trace_counter,:) =  ilxl_read(start_idx,:);
        results_out{1,2}{2,3}(trace_counter,:) =  uint32(1);     
         
        data = traces(:,start_idx:end_idx);
        trim_data = data;
        % Work out which way round the gathers are
        %if sum(trim_data(:,1) > 0) > sum(trim_data(:,end) > 0)
            trim_data_filt = fliplr(medfilt1(traces(:,start_idx:end_idx),5,[],2));
        %else
        %    trim_data_filt = medfilt1(traces(:,start_idx:end_idx),5,[],2);
        %end   
   
        [~,n_traces_gather] = size(trim_data);

        time = repmat((0:sample_rate:(n_samples-1)*sample_rate)',1,n_traces_gather);
        f_max = (1/sample_rate)*1000;

        S = (1/zsmooth)*spdiags(repmat([(1:1:zsmooth),(zsmooth-1:-1:1)],n_samples,1),[(-zsmooth+1:1:0),(1:1:zsmooth-1)],n_samples,n_samples);

        for ii = 2:n_traces_gather

            count = 0;

            t1 = trim_data_filt(:,ii-1);
            T1 = S*spdiags(t1,0,n_samples,n_samples);

            t2 = trim_data_filt(:,ii);

            phase_shifts = bsxfun(@times,(-shift:1:shift),(1/1000).*2.*pi.*repmat((0:f_max/(n_samples-1):f_max)',1,1+2*shift));
            t2_shifts = ifft(bsxfun(@times,fft(t2),exp(1i*phase_shifts)),'symmetric');

            for kk = -shift:1:shift

                count = count+1;

                T2 = S*spdiags(t2_shifts(:,count),0,n_samples,n_samples);

                det_coef(:,count) =  ((T1*t2_shifts(:,count)).*(T2*t1))./((T1*t1).*(T2*t2_shifts(:,count)));     

            end

            [~,idx] = max(det_coef');
            trim_shift(:,ii-1) = idx'-shift-1;

        end

        trim_shift = [zeros(n_samples,1),cumsum(fliplr(trim_shift),2)];
        trim_shift(data==0) = 0;
        
%         trim_sum = sum(trim_shift,2);
        n = length(trim_shift);
        trim_sum = sqrt((1/n)*sum((trim_shift.^2),2));
        
        for ii = 1:n_traces_gather
            trim_data(:,ii) = interp1(time(:,ii),data(:,ii),time(:,ii)-trim_shift(:,ii),'linear',0);
        end
        
        results_out{2,2}(:,start_idx:end_idx) = trim_shift;
        results_out{3,2}(:,start_idx:end_idx) = trim_data;
        results_out{4,2}(:,trace_counter) = trim_sum;
        
        trace_counter = trace_counter + 1;
        fprintf('%d complete\n',trace_counter);
        start_idx = start_idx + n_traces_gather;
        end_idx = start_idx;
        clear trim_sum trim_shift trim_data trim_data_filt
        end
        
    end   
    
    sw_segy_write_lite_gathers(results_out,i_block,sample_rate,output_dir);   
    
    %remove_job_directory(job_dir);    
      
end