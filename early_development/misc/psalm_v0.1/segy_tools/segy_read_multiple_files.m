function [traces positions] = segy_read_multiple_files(block_mat_all,process_files_mat,process_positions)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    process_files = load(process_files_mat,'-mat');

    for ii=1:process_files.nfiles
        [traces{ii} positions] = node_segy_read_traces(block_mat_all,process_positions,process_files_mat,ii);
    end

end

