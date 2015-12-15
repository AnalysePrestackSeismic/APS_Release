function input_files = njobs_worker(input_files,max_workers,n_blocks,n_offsets)
    
    % this function works with multiple files
    n_files = length(input_files);
    
    for i_files=1:1:n_files
        % have to make sure that we do not mix traces within aperture
        %aperture = input_files{i_files}.aperture;
        n_pos = n_offsets;%(2*aperture+1)^2;

        % traces in input volume    
        ntraces = input_files{i_files}.n_traces/n_offsets;

        % number of traces each worker will process
        ntraces_per_worker = floor(ntraces/max_workers);

        % number of traces each block will process   
        ntraces_per_block = floor(ntraces_per_worker/n_blocks);

        % number of traces left over - added to last block
        ntraces_left_over = ntraces - (ntraces_per_block*max_workers*n_blocks);

        % job index number
        job_i = 1:1:max_workers*n_blocks;

        % repeat matrix for total number of jobs; this is a vector
        ntraces_in_block = repmat(ntraces_per_block,1,n_blocks*max_workers);
        % add back remaining traces to last block
        ntraces_in_block(end) = ntraces_in_block(end)+ntraces_left_over;

        % calculate end index for each block
        end_index = job_i.*ntraces_per_block;         
        end_index(end) = end_index(end)+ntraces_left_over;

        % block index number
        block_i = 1:1:n_blocks;
        block_i = repmat(block_i',1,max_workers);

        % woker index number
        worker_i = 1:1:max_workers;
        worker_i = repmat(worker_i,n_blocks,1);

        % take account of aperture
        end_index = end_index.*n_pos;
        start_index = end_index-(ntraces_in_block.*n_pos)+1;

        %  outputs from function
        input_files{i_files}.index_worker = [job_i' worker_i(:) block_i(:) start_index' end_index'];    
        input_files{i_files}.job_info = [ntraces ntraces_per_worker ntraces_per_block ntraces_left_over]; 
    end
    
end