function RMSE_values = sensitivity_real_data(params)
% Sensitivity analysis by comparing with real data

    close all;
    reset_path;

    Dataset = parquetread('Data/Examples/Raj2020_Tests.parquet');

    % % Define the baseline parameters
    % load('params_EHM.mat');
    baseline_params = params;
    baseline_params.verbose = false;
    
    % Run the baseline simulation
    ModelName = 'EHM';
    Estimator = 'PEM';
    Target = 'Compare';

    addpath(genpath(strcat('./Code/Models/', ModelName)));
    addpath(genpath(strcat('./Code/Methods/', Estimator)));
    
    [Model, baseline_params] = step0(ModelName, 3, baseline_params);
    baseline_sol = step1(Target, Model, baseline_params, 3, Dataset);
    
    % save('baseline_voltage.mat', 'baseline_voltage');
    
    % Define parameters to vary and their range
    param_names_to_vary = {'In', 'Ip', 'Q', 'Rf', 'tau', 'nu', 'miu', 'alph', 'b'};

    % Define the range of delta
    min_value = -0.15;
    max_value = 0.15;
    num_values = 50;
    delta = linspace(min_value, max_value, num_values);
        
    % Initialize RMSE results array
    RMSE_values = zeros(length(param_names_to_vary), length(delta));
    
    % Loop over each parameter
    for i = 1:length(param_names_to_vary)
        for j = 1:length(delta)
            % Modify the parameter value
            varied_params = baseline_params;
            varied_params.(param_names_to_vary{i}) = baseline_params.(param_names_to_vary{i}) * (1 + delta(j));
            
            params.verbose = false;
            params.plot_results = false;

            % Run simulation with varied parameters
            [Model, params] = step0(ModelName, 3, varied_params);
            % pred_sol = run_simulation(Model, params);
            % Load or generate data

            [true_sol, params] = step1('Simulate', Model, params, 3, Dataset);
            
            % Perform estimation and update parameter values
            %[est_sol,  params] = step2(Target, Model, params, 3);
            
            % Run simulation using updated parameters
            %[pred_sol, params] = step3(Target, Model, params, 3, est_sol);
 
            params = step4('Compare', params, baseline_sol, true_sol);
            % Compute RMSE
            % Compute difference between data and simulation
            % Unpack RMSE Values
            RMSE_values(i, j) = params.RMSE_mV;
        end
    end
    
    figure; hold on;
    for i = 1:length(param_names_to_vary)
    plot(delta, RMSE_values(i, :), 'DisplayName', param_names_to_vary{i});
    end

    xlabel('Parameter Variation Delta');
    ylabel('RMSE (mV)');
    title('Sensitivity Analysis for EHM Parameters');
    legend('show');  % Displaying the legend
    hold off;

    % Display results
    disp(sensitivity_table);
end
