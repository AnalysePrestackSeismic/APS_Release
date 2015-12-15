function [seismic]=segy_trace_header_edit(inputfile,edit_byte,byte_value)
    
    tic

    seismic.filepath = inputfile;

    seismic.fid = fopen(seismic.filepath,'r+','ieee-be');

    fprintf('\nWriting %d to byte %d in trace headers of %s\n',byte_value,edit_byte,inputfile);

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
    
    % Overwrite first trace header byte in file
    fwrite(seismic.fid,byte_value,'int32',edit_byte-1);
    
    % Make array of new byte values
    byte_array = byte_value+zeros(seismic.n_traces-1,1);
    
    % Overwrite array of byte values in file
    fwrite(seismic.fid,byte_array,'int32',236+skip_trace);
    
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
