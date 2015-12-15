function [ ] = meta_data_2d_smo_xy( job_meta_path )
%ILXL_GRID Summary of this function goes here
% Input:
% Output:
%   Detailed explanation goes here
%%
job_meta = load(job_meta_path);

% minil = min(job_meta.block_keys(:,1));
% minxl = min(job_meta.block_keys(:,3));
% maxil = max(job_meta.block_keys(:,2));
% maxxl = max(job_meta.block_keys(:,4));

ils = sort(unique(job_meta.block_keys(:,1)));   % Sort inlines of  fist corner point of blocks in increasing order and remove duplicates
xls = sort(unique(job_meta.block_keys(:,3)));   % Sort xlines of  fist corner point of blocks in increasing order and remove duplicates

ilslen = length(ils);                           % Number of unique inlines of  fist corner point of blocks
xlslen = length(xls);                           % Number of unique inlines of  fist corner point of blocks

ilxl_lookup = zeros(ilslen,xlslen,2);           % Create look up table for inlines and cross-lines

loopfin = size(job_meta.liveblocks,1);          % Number of live blocks
lpi = 1;                                        % Loop Index
count = 1;                                      % Counter
%%
%Loop for creating look up table from Job Meta Information
while lpi <= loopfin
    i_block = job_meta.liveblocks(lpi);                                 % Block Number for Current Live Block
    if(i_block == 387)
        dave = 1;
    end
    currow = find(ils == job_meta.block_keys(i_block,1),1,'first');     % Find ROW = The inline matching ith block first coner point inline
    curcol = find(xls == job_meta.block_keys(i_block,3),1,'first');     % Find COULMN = The inline matching ith block first coner point inline
    
    ilxl_lookup(currow,curcol,1) = i_block;                             % Store Block Number in the first cell indexed by found ROW and COLUMN
    ilxl_lookup(currow,curcol,2) = job_meta.stdev(i_block);             % Store Block wavelet variance in the 2nd cell indexed by found ROW and COLUMN
    ilxl_lookup(currow,curcol,3) = job_meta.wb_z_avg(i_block);          % Store Block Average Water Bottom in the 3rd cell indexed by found ROW and COLUMN
    ilxl_lookup(currow,curcol,4) = job_meta.live_offset(i_block);       % Store Block Live Offset Value  in the 4th cell indexed by found ROW and COLUMN
    lpi = lpi + 1;                                                      % Increment Loop Index
end
%%
%Plotting
% check for NaNs
log_nan = isnan(ilxl_lookup(:,:,2));
if sum(sum(log_nan)) > 0
   ilxl_lookup_lin = squeeze(ilxl_lookup(:,:,2)); 
   ilxl_lookup_lin = ilxl_lookup_lin(:); 
   ilxl_lookup_lin(log_nan) = 0;
   ilxl_lookup(:,:,2) = reshape(ilxl_lookup_lin,xlslen,ilslen);   
end

figure(1); imagesc(ilxl_lookup(:,:,2)); title('input variance');%caxis([500 1000]);                     % Plot Variance
figure(11); imagesc(ilxl_lookup(:,:,3)); title('input water bottom depth'); %caxis([500 1000]);         % Plot Average Water Bottom Z
figure(12); imagesc(ilxl_lookup(:,:,4)); title('input live offset'); %caxis([0 50]);                    % Plot Live Offset
figure(14); imagesc(ilxl_lookup(:,:,1)); title('block number'); %caxis([0 50]);                         % Plot Block Number
%%
%Filtering
flen = 7;                                                                                           % Half Length of Filter
filt_smo2 =  [  linspace(0,flen,(flen+1)) , linspace((flen-1),0,flen)];                             % Smoothening Filter
filt_smo2 = filt_smo2/sum(filt_smo2);                                                               % Normalize Smoothening Filter
filt_smo = filt_smo2;

input_pad = ilxl_lookup(:,:,2);                                                                     % Initialize Input Pad to Variance Matrix ( il vs xl )

% Pad Inline and xline direction by the boundary values by (length of
% filter+1) cells
input_pad = [repmat(input_pad(:,1),1,(flen+1)), input_pad, repmat(input_pad(:,xlslen),1,(flen+1))]; 
input_pad = [repmat(input_pad(1,:),(flen+1),1); input_pad; repmat(input_pad(ilslen,:),(flen+1),1)];

ilxl_lookup_smo = conv2(filt_smo2,filt_smo,input_pad,'same');                                       % Smoothening : Convolve the designed filter with the padded input
ilxl_lookup_smoout = ilxl_lookup_smo((flen+2):(flen+2)+ilslen-1,(flen+2):(flen+2)+xlslen-1);        % Crop the smoothed output to original boundary
figure(2); imagesc(ilxl_lookup_smoout); title('smoothed variance');

ilxl_lookup_smoout(:,:,2) = ilxl_lookup(:,:,1);                             % Put Block Numbers the smoothed output as a second cell for every il xl
%figure(20); imagesc(ilxl_lookup_smoout(:,:,2)); title('block number');

ilxl_lookup_smoout = reshape(ilxl_lookup_smoout,[(ilslen*xlslen) 2]);       % Reshape into matrix indexed by Block Number, 1st cell: smoothed output, 2nd cell: bllock number for libe blocks

%ilxl_lookup_smo = sortrows(ilxl_lookup_smoout,2);

%figure(15); imagesc(ilxl_lookup_smo(:,:,1)); %caxis([0 50]);

job_meta.stdev_smo = ones(job_meta.n_blocks,1);                             % Initialized "to be" smoothed variance in job meta structure

% Loop for entering smoothed variance for live blocks
for ii = 1:size(ilxl_lookup_smoout)
    if ilxl_lookup_smoout(ii,2) > 0
        job_meta.stdev_smo(ilxl_lookup_smoout(ii,2)) = ilxl_lookup_smoout(ii,1);
    end
end

job_meta.live_offset_avg = round(mean(ilxl_lookup(currow,curcol,4)));

save(job_meta_path,'-struct','job_meta','-v7.3');


% ilxl_loc(:,1) = job_meta.block_keys(:,1) +  (((job_meta.block_keys(:,2) -  job_meta.block_keys(:,1))/job_meta.pkey_inc)/2)*job_meta.pkey_inc;
% ilxl_loc(:,2) = job_meta.block_keys(:,3) +  (((job_meta.block_keys(:,4) -  job_meta.block_keys(:,3))/job_meta.skey_inc)/2)*job_meta.skey_inc;
end