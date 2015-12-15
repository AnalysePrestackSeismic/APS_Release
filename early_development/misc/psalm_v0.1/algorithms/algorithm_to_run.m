function [distribute_type varargout] = algorithm_to_run(varargin)

    if nargin == 0
        fprintf('\nBased on your input files the following\n algorithms are available:\n');

        % Descriptive name, function name (path can found by using which)
        algorithms = {'View section' 'view_section';...
            'View slice' 'view_slice';...
            'RMS velocity inversion' 'rms_velocity_inversion';...   
            'Water Bottom Flatten' 'water_bottom_flatten';...            
            'Intercept Gradient Inversion' 'int_grad_inv_proj';...
            'Seismic anomaly spotter' 'seismic_anomaly_spotter';...
            'Automatic Fault Connector' 'automatic_fault_connector';...
            'Anomalous Body Connector' 'anomalous_body_connector';...
            'Wavelet Estimation' 'wavelet_estimation';...
            'Minimum Energy Chi' 'minimum_energy_chi';...
            
            };
            
            
%             'Constrained dix inversion' 'vrms_vint_inversion';...
%             'Wavelet estimation' 'wavelet_estimation';...
%             'Intercept gradient inversion' 'int_grad_inversion';...
%             'Simultaneous prestack inversion' 'sim_pre_stack_inversion';...
%             'Constrained water column fitter' 'water_column_fitter';...
%             'Phase estimation' 'bliss_kpe';...
%             'MCMC deconvolution' 'mcmc_deconvolution'
            

        for i_no = 1:1:length(algorithms)
           fprintf('[%d] %s \n',i_no,algorithms{i_no});
        end

        index = input('Enter the number of algorithm to run: ');
        distribute_type = algorithm_to_run(algorithms{index,2});
        varargout = {algorithms{index,1} algorithms{index,2}};

    elseif nargin > 0
        switch varargin{1}
            case 'view_section'
                distribute_type = 'N/A';
                % Ask section number
                
            case 'view_slice'
                distribute_type = 'N/A';
                % ask slice number
                
            case 'water_bottom_flatten' 
                distribute_type = 'trace';
                varargout = {'required_segy_types','Full stack','Angle stack'}; 
                
            case 'rms_velocity_inversion' 
                distribute_type = 'trace';
                varargout = {'required_segy_types','Full stack','Angle stack'};     
            
            case 'int_grad_inv_proj'
                distribute_type = 'trace';
                varargout = {'required_segy_types','Angle stack'}; 
                
            case 'wavelet_estimation'
                distribute_type = 'either';
                varargout = {'required_segy_types','Full stack','Angle stack'}; 
                
            case 'minimum_energy_chi'
                distribute_type = 'trace';
                varargout = {'required_segy_types','Angle stack'}; 
                
            case 'seismic_anomaly_spotter' 
                distribute_type = 'either';
                varargout = {'required_segy_types','Other attribute'...
                'required_algorithm','int_grad_inv_proj'...
                'input_flatten','Yes'
                }; 
                      
            case 'anomalous_body_connector' 
                distribute_type = 'slice';
                varargout = {'required_segy_types','Other attribute' 
                'required_algorithm','seismic_anomaly_finder'...
                };
            
            case 'automatic_fault_connector' 
                distribute_type = 'trace';
                varargout = {'required_segy_types','Other attribute' 
                'required_algorithm','seismic_anomaly_finder'...
                };
        end
    end

    % could have switch


end