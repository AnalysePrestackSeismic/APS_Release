function [slices] = optslice_segy(seis)
% Uses the output structure of segy_read_files to reorder the scanned segy
% file into fast slice access order. The output binary file is float32
% format.
% Jonathan Edgar 15/08/2012

for ii = 1:length(seis)

    % read and write in data blocks of 1GB
    block = floor((1024^3)/(4*seis{ii}.n_samples));
    loop = floor(seis{ii}.n_traces/block);

    % make an output file of zeros to be overwritten by the output data
    num = double(seis{ii}.filepaths);
    num = num == 47;
    num = cumsum(fliplr(num));
    num = num == 0;
    num = fliplr(num);
    fid_out = strcat(seis{ii}.filepaths(num),'_slices.bin');
    fid = fopen(fid_out,'a');
    for j=1:seis{ii}.n_samples
        fwrite(fid,zeros(seis{ii}.n_traces,1),'float32');
    end
    fclose(fid);

    % open the output file for overwriting
    fid = fopen(fid_out,'r+');

    % loop to read blocks of data traces and write out in slice order
    for i=1:loop
        fprintf('Reading traces %d to %d of %d\n',1+(i-1)*block,i*block,seis{ii}.n_traces);
        [traces] = segy_read_traces(seis{ii},1+(i-1)*block,i*block,0,0);
        fseek(fid,(i-1)*block*4,'bof');

        fprintf('Writing traces %d to %d of %d\n',1+(i-1)*block,i*block,seis{ii}.n_traces);
        for k=1:seis{ii}.n_samples
            fwrite(fid,traces.data(k,:),'float32');
            fseek(fid,(seis{ii}.n_traces-block)*4,'cof');
        end
    end

    % work out the number of traces yet to read
    leftovers = seis{ii}.n_traces-(loop*block);

    % read the remaining block of data traces and write out in slice order
    fprintf('Reading traces %d to %d of %d\n',1+(loop*block),seis{ii}.n_traces,seis{ii}.n_traces);
    [traces] = segy_read_traces(seis{ii},1+(loop*block),seis{ii}.n_traces,0,0);
    fseek(fid,loop*block*4,'bof');

    fprintf('Writing traces %d to %d of %d\n',1+(loop*block),seis{ii}.n_traces,seis{ii}.n_traces);
    for k=1:seis{ii}.n_samples
        fwrite(fid,traces.data(k,:),'float32');
        fseek(fid,(seis{ii}.n_traces-leftovers)*4,'cof');
    end

    fclose(fid);

    slices{ii}.n_traces = seis{ii}.n_traces;
    slices{ii}.n_samples = seis{ii}.n_samples;
    slices{ii}.min_iline = seis{ii}.min_iline;
    slices{ii}.max_iline = seis{ii}.max_iline;
    slices{ii}.min_xline = seis{ii}.min_xline;
    slices{ii}.max_xline = seis{ii}.max_xline;
    slices{ii}.n_iline = seis{ii}.n_iline;
    slices{ii}.n_xline = seis{ii}.n_xline;
    slices{ii}.il_inc = seis{ii}.il_inc;
    slices{ii}.xl_inc = seis{ii}.xl_inc;
    slices{ii}.slice_pointers = [(1:1:slices{ii}.n_samples)',(0:slices{ii}.n_traces*4:slices{ii}.n_traces*(slices{ii}.n_samples-1)*4)'];
    slices{ii}.geometry = [0;seis{ii}.s_rate/1e3;seis{ii}.n_samples;reshape(seis{ii}.trace_pointers(:,1:2)',2*seis{ii}.n_traces,[])];
    if isfield(seis{ii},'angle')
        slices{ii}.angle = seis{ii}.angle;
    end

end