function [th] = scansegy(file_in,varargin)

    % Open the segy file to read
    fid = fopen(file_in,'r','b');
    
    % Calculate the number of trace header bytes to read
    nbytes_to_read = size(varargin,2);
    
    % Allocate memory for indices of bytes to read
    idx_to_read = zeros(nbytes_to_read,2);
    
    % Calculate the indicies of the trace header bytes to read and flag if they need co-ordinate scaler adjustment
    for ii = 1:nbytes_to_read
        idx_to_read(ii,1) = 1+((varargin{ii}-1)/4);
        if or(and(varargin{ii}>=73,varargin{ii}<=85),and(varargin{ii}>=181,varargin{ii}<=185))
            idx_to_read(ii,2) = 1;
        end
    end

    % Skip text header
    fseek(fid,3200,'bof');

    % Read binary header as 
    bh=fread(fid,400,'uint8');

    % Re-interpret binary header as uint16 or uint32 as required
    two_bytes=bh(1:2:399)*256+bh(2:2:400);
    four_bytes=((bh(1:4:9)*256+bh(2:4:10))*256+bh(3:4:11))*256+bh(4:4:12);
    bh=[four_bytes(1:3);two_bytes(7:200)];

    % Get number of samples per trace
    nsamples = bh(8);
    
    % Calculate number of traces
    ll=dir(file_in);
    ntr=0.25*(ll.bytes-3600)/(nsamples+60);
    
    % Allocate memory for trace headers to be read
    th = zeros(nbytes_to_read,ntr,'uint32');
    
    % Set number of traces per block (reduce for low memory systems)
    blocktr = 1000;
    loop_end = floor(ntr/blocktr);

    % Loop to read segy data
    tic
    for ii = 1:loop_end
        % Read blocktr x trace headers and trace data as uint32 into a temporary matrix
        tmptr = fread(fid,[60+nsamples,blocktr],'uint32=>uint32');
        if ii == 1
            % Re-interpret bytes 71-72 as uint16 and convert from two's compliment
            coscaler = bin2dec(num2str(bitget(tmptr(18,1),(16:-1:1))))-2^16;
        end
        % Extract only the trace headers requested by varargin
        th(1:nbytes_to_read,1+(ii-1)*blocktr:ii*blocktr) = tmptr(idx_to_read(:,1),:);
    end
    
    % Calculate the number of trace headers not read by the loop above
    leftovers = ntr-loop_end*blocktr;
    
    % Read the remaining trace headers (if any)
    if leftovers > 0
        tmptr = fread(fid,[60+nsamples,leftovers],'uint32=>uint32');
        th(1:nbytes_to_read,1+loop_end*blocktr:ntr) = tmptr(idx_to_read(:,1),:);
    end
    toc
    fprintf('%d MB of data read at %d MB/sec in %d seconds\n',round((ll.bytes-3600)/(1024*1024)),round(((ll.bytes-3600)/(1024*1024))/toc),round(toc));
    
    % Apply the co-ordinate scaler
    if coscaler < 0
        coscaler = 1/abs(coscaler);
        th(logical(idx_to_read(:,2)),:) = th(logical(idx_to_read(:,2)),:)*coscaler;
    elseif coscaler > 0
        th(logical(idx_to_read(:,2)),:) = th(logical(idx_to_read(:,2)),:)*coscaler;
    end
    
    fclose all;
end