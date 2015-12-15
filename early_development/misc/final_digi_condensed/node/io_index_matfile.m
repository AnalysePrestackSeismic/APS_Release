function [] = io_index_matfile()

    % Load seismic structre (scanned trace headers)
    % seismic = load(seismic_mat_path);
    % Check with full stack or gathers
       
    seismic_mat_path = '/TZA/segy/2013_kusini_inboard/pgs_standard_volume/angle_stacks/fordigi/05-10_angstk.mat_lite';
        tmp_seismic = fread(fopen(seismic_mat_path,'r'),'double');
        
        % read inline and crossline bytes from index file
        % read any other bytes too (e.g. vargargin for offset etc)
        seismic.filepath = char(tmp_seismic(1:2000,1)');
        seismic.file_type = tmp_seismic(2001,1);
        seismic.n_samples = tmp_seismic(2002,1);
        seismic.trace_ilxl_bytes = tmp_seismic(2003:end);
       
        
        outfile_mat = '/TZA/segy/2013_kusini_inboard/pgs_standard_volume/angle_stacks/fordigi/20-25_angstk.mat_lite';
        outfile_sgy = '/TZA/segy/2013_kusini_inboard/pgs_standard_volume/angle_stacks/fordigi/20-25_angstk.sgy';
        filepath_binary = uint64(outfile_sgy);
        pad_filepath = zeros(1,(2000-length(filepath_binary)));
        filepath_binary = [filepath_binary,pad_filepath];
        fwrite(fopen(outfile_mat,'w'),[filepath_binary';seismic.file_type;seismic.n_samples;seismic.trace_ilxl_bytes],'double');   
        
        
        outfile_mat = '/TZA/segy/2013_kusini_inboard/pgs_standard_volume/angle_stacks/fordigi/25-30_angstk.mat_lite';
        outfile_sgy = '/TZA/segy/2013_kusini_inboard/pgs_standard_volume/angle_stacks/fordigi/25-30_angstk.sgy';
        filepath_binary = uint64(outfile_sgy);
        pad_filepath = zeros(1,(2000-length(filepath_binary)));
        filepath_binary = [filepath_binary,pad_filepath];
        fwrite(fopen(outfile_mat,'w'),[filepath_binary';seismic.file_type;seismic.n_samples;seismic.trace_ilxl_bytes],'double');           
        
        
        outfile_mat = '/TZA/segy/2013_kusini_inboard/pgs_standard_volume/angle_stacks/fordigi/30-35_angstk.mat_lite';
        outfile_sgy = '/TZA/segy/2013_kusini_inboard/pgs_standard_volume/angle_stacks/fordigi/30-35_angstk.sgy';
        filepath_binary = uint64(outfile_sgy);
        pad_filepath = zeros(1,(2000-length(filepath_binary)));
        filepath_binary = [filepath_binary,pad_filepath];
        fwrite(fopen(outfile_mat,'w'),[filepath_binary';seismic.file_type;seismic.n_samples;seismic.trace_ilxl_bytes],'double');           
        
        outfile_mat = '/TZA/segy/2013_kusini_inboard/pgs_standard_volume/angle_stacks/fordigi/35-40_angstk.mat_lite';
        outfile_sgy = '/TZA/segy/2013_kusini_inboard/pgs_standard_volume/angle_stacks/fordigi/35-40_angstk.sgy';
        filepath_binary = uint64(outfile_sgy);
        pad_filepath = zeros(1,(2000-length(filepath_binary)));
        filepath_binary = [filepath_binary,pad_filepath];
        fwrite(fopen(outfile_mat,'w'),[filepath_binary';seismic.file_type;seismic.n_samples;seismic.trace_ilxl_bytes],'double');   
        
        
        outfile_mat = '/TZA/segy/2013_kusini_inboard/pgs_standard_volume/angle_stacks/fordigi/40-45_angstk.mat_lite';
        outfile_sgy = '/TZA/segy/2013_kusini_inboard/pgs_standard_volume/angle_stacks/fordigi/40-45_angstk.sgy';
        filepath_binary = uint64(outfile_sgy);
        pad_filepath = zeros(1,(2000-length(filepath_binary)));
        filepath_binary = [filepath_binary,pad_filepath];
        fwrite(fopen(outfile_mat,'w'),[filepath_binary';seismic.file_type;seismic.n_samples;seismic.trace_ilxl_bytes],'double');           
                
        
        fclose('all');  
end