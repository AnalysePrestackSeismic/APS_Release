% batch(sch,fit_2d_gauss,0,{input_segy{1}.filepaths,input_segy{1}.index_worker(start_batch_index:ii*n_blocks,:),n_freq,output_dir},'matlabpool',0) 

function fit_2d_gauss(filepath,n_samples,traces_process,n_freq,res_inc,output_location)
%final 

% Read traces for multiple gathers
traces = segy_read_traces(filepath,n_samples,traces_process);

% Fit 2d gauss to single gather
ntraces = size(traces.data,2);

start_trace = 1;
n_gathers = ntraces/n_freq;

f_info.start = 1;
f_info.end = n_freq;

t_info.start = 3000;
t_info.end = 3800;

res_min = 1;
num_mors = 30;

search.sigt=0.25;
search.sigw=0.25;


    for i=1:1:n_gathers

         if (traces.pos(1,1)==9106)
            vals=traces.data(:,start_trace:start_trace+n_freq-1);

        %xline_num=round((i-1)/n_freq+first_xl);
        %if xline==(first_xl+n_xlines); iline=iline+1; xline=first_xl;p=p+1;q=1;
        %elseif j==1;xline=first_xl;else xline=xline+1; end

        %%%%%scale values to between 0 and 255  
        scale=255/max(max(vals));
        vals_scale=vals*scale;

        %%%%%%upscale picture by a factor of res_inc
        vals2=iir_adapted(vals_scale,res_inc,'display','off');
        %%%%create new frequency and time axis for upscaled image
        freq2=linspace(f_info.start,f_info.end,size(vals2,2));
        time2=linspace(t_info.start,t_info.end,size(vals2,1));
        %%%% time and frequency sampling are altered by upscaling image
        t_info.sampling=(t_info.end-t_info.start)/size(vals2,1);
        f_info.sampling=(f_info.end-f_info.start)/size(vals2,2);

        param_store{i} = fit_2d_morlet_test2(time2,freq2,vals2,search,res_min,num_mors,traces.pos(1,start_trace),traces.pos(2,start_trace));

        else param_store{i}=1;
        end
        
        start_trace = start_trace + n_freq;
        

    end

        start_il = num2str(traces.pos(1,1));
        end_il = num2str(traces.pos(1,end));

        start_xl = num2str(traces.pos(2,1));
        end_xl = num2str(traces.pos(2,end));



    % Save results of fitting
        data_out = strcat(output_location,'/','param_store_',start_il,'-',end_il,'_xl',start_xl,'-',end_xl);
        save(data_out,'param_store');
        %file_vint = fopen(data_out,'a');
        %fwrite(file_vint,param_store,'float32');
        %fclose(file_vint);

        % Write position information
        % meta_out = strcat(output_location,'/','vrms_pos_worker_',worker,'_block_',block);
        % dlmwrite(meta_out,vint.pos','delimiter','\t');

        %figure(i_block)
        %imagesc(vint.m)
    
end