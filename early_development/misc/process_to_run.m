function [algo_name func_name] = process_to_run(print_all)

if print_all == 1
    fprintf('The following algorithms are available:\n');
    
    % Descriptive name, function name (path can found by using which)
    algorithms = {'Structural closure finder' 'path';...
        'Seismic anomaly hunter' 'path';...
        'Constrained dix inversion' 'vrms_vint_inversion';...
        'Intercept gradient inversion' 'int_grad_inversion';...
        'Simultaneous prestack inversion' 'sim_pre_stack_inversion';...
        'Constrained water column fitter' 'water_column_fitter';...
        'Phase estimation' 'bliss_kpe';...
        'MCMC deconvolution' 'mcmc_deconvolution'};
    
    for i_no = 1:1:length(algorithms)
       fprintf('[%d] %s \n',i_no,algorithms{i_no});
    end

    index = input('Enter the number of algorithm to run: ');
    algo_name = algorithms{index,1};
    func_name = algorithms{index,2};
else
    



end