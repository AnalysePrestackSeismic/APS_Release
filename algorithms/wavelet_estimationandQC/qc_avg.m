function [] =qc_avg(job_meta_path)
%% FUNCTION FOR AVERAGING PRE DIGI QC OVER ALL BLOCKS
% Plots AVA of RMS Amplitude,
% Plots AVA of Variance
%%
close all;
job_meta = load(job_meta_path);


loopfin = size(job_meta.liveblocks,1);      % Number of Live Blocks
lpi = 1;                                    % Loop Index
count=1;
count = 1;                                  % Counter for blocks that have a corresponding wavelet file
ava_qc_all = struct;
while lpi <= loopfin
    i_block = job_meta.liveblocks(lpi);     % Reference the block numbers for live blocks
    qc_filepath = strcat(job_meta.ava_qc_directory,'ava_qc_',num2str(i_block),'.mat');
    tmpfileout=sprintf(strcat('opening file : ',qc_filepath,'\n'));
    
    if exist(qc_filepath,'file')
        disp(tmpfileout);
        
        ava_qc=load(qc_filepath);
        ava_qc_all.i_vol_max(count)=ava_qc.i_vol_max;
        ava_qc_all.n_win(count)=ava_qc.n_win;
        ava_qc_all.pick_wb_ind(count)=ava_qc.pick_wb_ind;
        ava_qc_all.n_traces(count)=ava_qc.n_traces;
        ava_qc_all.wb_ind_avg(count)=ava_qc.wb_ind_avg;
        
        n_vol=ava_qc.i_vol_max;
        n_win=ava_qc.n_win;
        if count == 1                                           % Initialize the wavelet matrix if this is the first live block
            ava_qc_all.rms = zeros(n_vol,n_win,loopfin);        % Matrix for storing all wavelets for all volumes for all live blocks, (Note: length of wavlet array = number of samples in wavelet) + 2
            ava_qc_all.rms(:,:,count) = ava_qc.rms;             % Reshape the Matrix sperating the different angle volumes in different columns for the first block              
            max_n_win = n_win;                                  % Intitialize maximum number of wavelet windows in block as number of windows in first block
            max_n_vol = n_vol;                                  % Intitialize maximum number of angle volumes in block as number of windows in first block
            total_n_traces= ava_qc.n_traces;       
        else                                                    % If this is not the first live block, append already initialized wavelet matrix
            if n_win > max_n_win                                % If you encounter a block with bigger window than what has been encountered ljust append a slab for the extra window
                rms_append = zeros(n_vol,(n_win-max_n_win),loopfin); % Create new slab of the extra length of window for all volumes for all live blocks
                ava_qc_all.rms = [ava_qc_all.rms rms_append];                      % Append the new slab for sccomodating the extra length of window
                ava_qc_all.rms(:,1:n_win,count) = ava_qc.rms; % Reshape the Matrix sperating the different angle volumes in different columns for the current block 
                max_n_win = n_win;
            else                                                % If this block has the same or less number of wavelet windows than encountered before 
                ava_qc_all.rms(:,1:n_win,count) = ava_qc.rms;% Reshape the Matrix sperating the different angle volumes in different columns for the current block
            end
            total_n_traces=total_n_traces + ava_qc.n_traces;
        end
        lpi = lpi + 1;                                          % Increment loop index 
        count  = count + 1;                                     % Increment counter
        
    else                                                        % If file doesnot exist
        tmpfileout = sprintf('Error opening file %s',qc_filepath);
        disp(tmpfileout)                                        % Display error in opening file
        lpi = lpi + 1;                                          % Increment loop index 
    end
    
end
count = count -1;
sum=0;
for ii=1:count
    sum = sum+ava_qc_all.n_traces(ii)*(ava_qc_all.rms(:,:,ii).*ava_qc_all.rms(:,:,ii));
end

ava_qc_all.rms_rms= sqrt(sum)/total_n_traces;                                       % Sum the wavelet spectrum including number of live traces in the blocks for all the blocks

angles=job_meta.tkey_min:job_meta.tkey_inc:job_meta.tkey_max;
figure(1);
for k=1:max_n_win
    t=strcat('level: ',num2str(((k-0.5)*job_meta.ns_overlap_qc*job_meta.s_rate/1000)),'m/ms');
    subplot(1,max_n_win,k);
    scatter(angles,ava_qc_all.rms_rms(:,k),'MarkerEdgeColor','b','MarkerFaceColor','c');xlim([min(angles) max(angles)]);
    xlabel( 'angle'); ylabel( 'RMS Amplitude');
    title (t);
end

end


