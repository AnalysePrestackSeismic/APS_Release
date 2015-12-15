function [] = optil_sw(sliceoptfile_in,slice_location,seismic_mat_path,bytespersample,mem,i_block,n_blocks)
% Optimises binary files for inline access and recombines geometery information so that the output can be read by opendtect.
%
% USAGE:   [] = optil(sliceoptfile_in,geomfile_in,bytespersample,mem)
%
% INPUTS:   sliceoptfile_in = slice optimised binary file as output by the function optslice.
%           geomfile_in = geometry file as output by the function geomstrip.
%           bytespersample = number of bytes per sample (4 if input data is float32 format).
%           mem = is the memory block in Gb for use by the code. Actual memory used will be higher than this.
%
% OUTPUTS:  inlineopt_sliceoptfile_in = inline optimised binary file including sampling info and inline and crossline loacions.




%% To do: 
% Best to run the optil function on one file so built a system call in a
% similar way to the srun command to CAT all of the SAS binary files
% together. Make this so it writes out the full file list in the right
% order and then submits it as a system call to linux.








%% Srun code to adapt

for i_block = 1:1:n_blocks
    file_array{1,i_block} = sprintf('%ssas_result_slices_block_%d.bin',slice_location,i_block);
    
    %cat_command = sprintf('cat %ssas_result_slices_block_%d.bin',slice_location,i_block)
end

cat_this = sprintf('%s ',file_array{:});

cat_command = sprintf('cat %s > %ssas_result_volume.bin',cat_this,slice_location)

system(cat_command);


         slurm_job = sprintf('srun -c 2 -p %s -J %s_block_%d %s %s %s %s %s %s %s %s',slurm_part,algorithm_name,i_block,run_script_path,matlab_path,algorithm_name,seismic_mat_path,num2str(i_block),num2str(n_blocks),arguments,'&');
         %slurm_job = sprintf('srun -c 12 -p %s -J %s_block_%d %s %s %s %s %s %s %s %s',slurm_part,algorithm_name,i_block,run_script_path,matlab_path,algorithm_name,seismic_mat_path,num2str(i_block),num2str(n_blocks),arguments,'&');
    %    slurm_job = sprintf('srun -x /apps/gsc/matlab-mcode-beta/eslib/psalm_lite/node/nodelist -c 1 -p %s -J %s_block_%d %s %s %s %s %s %s %s %s',slurm_part,algorithm_name,i_block,run_script_path,matlab_path,algorithm_name,seismic_mat_path,num2str(i_block),num2str(n_blocks),arguments,'&');

        system(slurm_job);
        fprintf('Function %s, job %d submitted to SLURM\n',algorithm_name,i_block);
     


















%% Optil code.

% Open slice optimised input file.
fid1 = fopen(sliceoptfile_in);

% Open geometry file.
fid2 = fopen(geomfile_in);

% Create output inline optimised file containing geometry info.
fid3 = fopen(sprintf('%s_inlineopt',sliceoptfile_in),'a');

% Work out the size of the slice optimised input file.
fseek(fid1, 0, 'eof');
filesize = ftell(fid1);
frewind(fid1);

% Write start time to output file.
fwrite(fid3,fread(fid2,1,'*float32'),'float32');

% Write sample rate to output file.
fwrite(fid3,fread(fid2,1,'*float32'),'float32');

% Read number of time samples from geometry file.
nt = fread(fid2,1,'*int');

% Write number of time samples to output file.
fwrite(fid3,nt,'int');

% Change type from int to double to avoid divide errors.
nt = double(nt);

% Calculate the number of traces in the slice optimised input file.
ntrace = filesize/(nt*bytespersample);

% Work out the block size in time samples to read such that the size of the block is equal to 'mem'.
block = mem*(1024^3)/(nt*bytespersample);

% Ensure the block size is smaller than the number of time samples.
if block >= ntrace
    block = ntrace;
end

% Initialise some variables.
loops = floor(ntrace/block);
block = floor(block);
tmpil = zeros(nt,block);

% Loops to read and write the data.
for j = 1:loops
    % Print read progress to screen.
    fprintf('reading traces %d to %d of %d\n',1+((j-1)*block),j*block,ntrace)
    for k = 1:nt
        % Read block of traces at constant time slice.
        tmpslice = fread(fid1,block,'*float32');
        % Jump to same block of traces on the next time slice down.
        fseek(fid1,(ntrace-block)*bytespersample,'cof');
        % Build matrix of traces read.
        tmpil(k,:) = tmpslice';
    end
    % Print write progress to screen.
    fprintf('writing traces %d to %d of %d\n',1+((j-1)*block),j*block,ntrace)
    % Loop over all traces in the block.
    for k = 1:block
        % Write trace geometry in block to file.
%         fwrite(fid3,fread(fid2,2,'*int'),'int');
        % Write traces in block to file.
        fwrite(fid3,tmpil(:,k),'float32');
    end
    % Jump to next block of traces to read.
    fseek(fid1,(j*block)*bytespersample,'bof');
end

% Loops to read the remaining part of the data.        
if block*loops ~= ntrace
    leftovers = ntrace - (block*loops);
    tmpil = zeros(nt,leftovers);
    fprintf('reading traces %d to %d of %d\n',j*block+1,ntrace,ntrace)
    for k = 1:nt
        tmpslice = fread(fid1,leftovers,'*float32');
        fseek(fid1,(ntrace-leftovers)*bytespersample,'cof');
        tmpil(k,:) = tmpslice';
    end
    fprintf('writing traces %d to %d of %d\n',j*block+1,ntrace,ntrace)
    for k = 1:leftovers
%         fwrite(fid3,fread(fid2,2,'*int'),'int');
        fwrite(fid3,tmpil(:,k),'float32');
    end
end

% Close all files.
fclose all;

end
