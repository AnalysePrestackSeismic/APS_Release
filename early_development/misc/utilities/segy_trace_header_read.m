function [seismic]=segy_trace_header_read(inputfile,read_byte)
    
    tic

    seismic.filepath = inputfile;

    seismic.fid = fopen(seismic.filepath,'r','ieee-be');

    fprintf('\nReading %d byte in trace headers of %s\n',read_byte,inputfile);

    % Skip textual header (3200 bytes)
    fseek(seismic.fid,3200,'bof');

    % Read binary header (400 bytes)
    seismic.binary_header(:,:) = read_binary_file_header(seismic.fid);

    % Enter number of samples from binary header into seismic structure
    seismic.n_samples = seismic.binary_header(8);

    % Find out number of traces
    seismic.n_traces = get_no_of_traces(seismic.filepath,seismic.n_samples);

    % Calculate number of trace data bytes to skip
    skip_trace = seismic.n_samples*4;
    
    % Skip to first trace header byte to read in file
    fseek(seismic.fid,read_byte-1,'cof');
        
    % Read byte values in file into an array
    seismic.byte_array = fread(seismic.fid,[seismic.n_traces,1],'int32=>int32',236+skip_trace);
    
    toc
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function binary_header=read_binary_file_header(fid)
% Read binary header and initiate the seismic dataset

    bh=fread(fid,400,'uchar');

    two_bytes=bh(1:2:399)*256+bh(2:2:400);
    four_bytes=((bh(1:4:9)*256+bh(2:4:10))*256+bh(3:4:11))*256+bh(4:4:12);

    binary_header=[four_bytes(1:3);two_bytes(7:200)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ntr=get_no_of_traces(filename,nsamples)
% Use number of bytes to compute number of traces

    ll=dir(filename);
    nbytes=ll.bytes;
    ntr=0.25*(nbytes-3600)/(nsamples+60);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
