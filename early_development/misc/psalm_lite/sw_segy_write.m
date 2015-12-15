%function [] segy_write(input_dir, filename)


%Combine ILXL information with trace data <----add this to end of inversion
%code


%[(results_out{4,2}{1,2})';results_out{1,2}]



start_time = 0;
sample_rate = 4;
coscal_val = -100;
block_sz = 100; % number of traces written at a time

n_samples = size(results_out{1,2},1);

% open the file for writting
fid_ilxl_f32 = fopen('sw_test12.sgy','w');
%fid_ilxl_f32 = fopen('sw_test8.sgy','w','ieee-le');

header = zeros(60,block_sz,'int32');
%arrayofones = ones(1,block_sz,'int32');
% two lots of 16ints first set to 1 written as a 32, second set to zero.
% i.e when read as 16bit will come out as 1 0 
%arrayofones65 = arrayofones*65536;


%set an array for the coscal as packed int16 into int32 data
tmpcolscal = zeros(1,(block_sz*2),'int16');
tmpcolscal(1:2:size(tmpcolscal,2)) = coscal_val;
tmpcolscal2 = typecast(tmpcolscal,'int32'); 

%Define binary header
bin_header = zeros(200,1,'uint16');
bin_header(7) = 1;     % number of data traces per ensemble
bin_header(9) = sample_rate*1000; %Sample rate needs to be in microseconds
bin_header(11) = n_samples;  % number of samples per trace
bin_header(13) = 5;      % segy data number format 5 = ieee float 754
bin_header(14) = 1;      % ensemble fold
bin_header(15) = 4;      % trace sorting code 4 = horizontally stacked , 2 = cdp ensemble
bin_header(28) = 1;     % measurement system
bin_header(151) = 256;   %segy revision number
bin_header(152) = 1;      % fixed length flag  1= size always the same 


% write out the ebcdic_header
seismic.text_header = ['AAAAA',blanks(3190),'ZZZZZ'];
%seismic.text_header = char(ebcdic2ascii(reshape(seismic.text_header,80,40)')); % REVERSE THIS
fwrite(fid_ilxl_f32,seismic.text_header,'char*1',0,'ieee-be');

% write out the binary header

%     %GENERATE 400BYTE STRING TO POPULATE
%     % Re-interpret binary header as uint16 or uint32 as required
%     two_bytes=seismic.binary_header(1:2:399)*256+seismic.binary_header(2:2:400);
%     four_bytes=((seismic.binary_header(1:4:9)*256+seismic.binary_header(2:4:10))*256+seismic.binary_header(3:4:11))*256+seismic.binary_header(4:4:12);
%     seismic.binary_header=[four_bytes(1:3);two_bytes(7:200)];
%
fwrite(fid_ilxl_f32,bin_header,'uint16',0,'ieee-be');

% now loop through the data and write out

%count = [1:1:block_sz];


tic

for ii= 1:block_sz:size(results_out{1,2},2)
    
    maxblock = ii + block_sz - 1;
    
    if maxblock >  size(results_out{1,2},2)
        maxblock = size(results_out{1,2},2);
        block_sz = (maxblock - ii+1);
        header = zeros(60,block_sz,'int32');
       % arrayofones = ones(1,block_sz,'int32');
       % arrayofones65 = arrayofones*65536;
        % set the 
        %cjtmp = typecast(arrayofones*coscal_val,'int16');
        %cjtmp(2:2:size(cjtmp,2)) = 0;
        %cjtmp2 = typecast(cjtmp,'int32'); 
        %set an array for the coscal as packed int16 into int32 data
        tmpcolscal = zeros(1,(block_sz*2),'int16');
        tmpcolscal(1:2:size(tmpcolscal,2)) = coscal_val;
        tmpcolscal2 = typecast(tmpcolscal,'int32'); 
    end
    
    %temparr = single(results_out{1,2}(:,ii:maxblock));
    
    
    
    %for n = 1:block_sz
    
    %temparr(:,:) = typecast(swapbytes(single(results_out{1,2}(:,ii:maxblock))),'int32');
    
    %%temparr = typecast(swapbytes(single(reshape(results_out{1,2}(:,ii:maxblock),1,(n_samples*(maxblock-ii+1))))),'int32');
    temparr = typecast(single(reshape(results_out{1,2}(:,ii:maxblock),1,(n_samples*(maxblock-ii+1)))),'int32');
    %temparr = typecast(swapbytes(single(reshape(results_out{1,2}(:,ii:maxblock),1,(n_samples*(maxblock-ii+1))))),'int32');
    
    %temparr2 = reshape(temparr,n_samples,block_sz);
    
    %temparr2 = [header; swapbytes(reshape(temparr,n_samples,block_sz))];
    temparr2 = [header; reshape(temparr,n_samples,block_sz)];
    %need to resize last block as all are 100 long at the moment.
    
    %Set header values
    count = ii:1:maxblock;
    temparr2(1,:) = count;
    
    % set byte 29 to 1 as a 16 bit integer into a int32, so the first 16
    % bytes have to represent 1 when read as unit16, so multiply 1 by 2^16
    % which is 65536
    temparr2(8,:) = 65536; %arrayofones65; 
    
    % to write the 16int coscal as -100 as int32 in to byte location 71
    temparr2(18,:) = tmpcolscal2;

    
    %temparr2(47:48,:) = (results_out{4,2}{1,1}(ii:maxblock,:)';
    % write the inline crossline numbers to bytes 189 and 193
    temparr2(48,:) = typecast(results_out{4,2}{1,1}(ii:maxblock,1)','int32');
    temparr2(49,:) = typecast(results_out{4,2}{1,1}(ii:maxblock,2)','int32');
    
    %end
    
    
    fwrite(fid_ilxl_f32,temparr2,'int32',0,'ieee-be');
    
    %
    %
    %
    %
end
fclose(fid_ilxl_f32);

'Time to write file' = toc