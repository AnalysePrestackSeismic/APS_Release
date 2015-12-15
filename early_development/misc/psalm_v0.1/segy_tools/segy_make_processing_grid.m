function processing_grid = segy_make_processing_grid(ilxl_step,ilxl_aperture,ilxl_aperture_step,z_step,z_aperture,process_files)
           
    % Work out grid steps
    n_il_steps = ilxl_step*process_files.il_inc;
    n_xl_steps = ilxl_step*process_files.xl_inc;

    % set up grid
    grid_il = process_files.min_iline:n_il_steps:process_files.max_iline;
    grid_xl = process_files.min_xline:n_xl_steps:process_files.max_xline;

    % replicate matrices for inline and crosslines positions
    grid_xl_rep = repmat(grid_xl,1,length(grid_il));
    grid_il_rep = sort(repmat(grid_il,1,length(grid_xl)));

    % create z positions grid
    z_positions_up = floor(process_files.n_samples/2):z_step:process_files.n_samples;
    z_positions_down = fliplr(floor(process_files.n_samples/2)-z_step:-z_step:1);

    processing_grid.z_grid = [z_positions_down'; z_positions_up'];
    processing_grid.ilxl_grid = [grid_il_rep' grid_xl_rep'];

    % save some variables
    processing_grid.ilxl_aperture = ilxl_aperture;
    processing_grid.z_aperture = z_aperture;
    processing_grid.z_step = z_step;
    processing_grid.ilxl_step = ilxl_step;
    processing_grid.ilxl_aperture_step = ilxl_aperture_step;
end