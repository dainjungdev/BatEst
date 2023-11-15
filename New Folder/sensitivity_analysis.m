function RMSE_values = sensitivity_analysis(params)
    close all;
    clear all;
    reset_path;

    Dataset = parquetread('Data/Examples/Raj2020_Tests.parquet');

    % Define the baseline parameters
    load('params_EHM.mat');
    baseline_params = params;
    baseline_params.verbose = false;
    
    % Run the baseline simulation
    ModelName = 'EHM';
    Estimator = 'PEM';
    Target = 'Simulate';

    addpath(genpath(strcat('./Code/Models/', ModelName)));
    addpath(genpath(strcat('./Code/Methods/', Estimator)));
    
    [Model, baseline_params] = step0(ModelName, 3, baseline_params);
    baseline_sol = step1(Target, Model, baseline_params, 3, Dataset);
    
    % save('baseline_voltage.mat', 'baseline_voltage');
    
    % Define parameters to vary and their range
    param_names_to_vary = {'In', 'Ip', 'Q', 'Rf', 'tau', 'alph'};
    % Q: not realistic value of voltage - warning?

    % Crate: safety limit & scale plots but don't change the simulation
    delta = [-0.2, -0.1, 0, 0.1, 0.2];
        
    % Initialize RMSE results array
    RMSE_values = zeros(length(param_names_to_vary), length(delta));
    
    % Loop over each parameter
    for i = 1:length(param_names_to_vary)
        for j = 1:length(delta)
            % Modify the parameter value
            varied_params = baseline_params;
            varied_params.(param_names_to_vary{i}) = baseline_params.(param_names_to_vary{i}) * (1 + delta(j));
            
            disp(baseline_params.Rf)
            disp(varied_params.Rf)


            % Run simulation with varied parameters
            [Model, params] = step0(ModelName, 3, varied_params);
            % pred_sol = run_simulation(Model, params);
            % Load or generate data
            params.verbose = false;
            params.plot_results = false;

            [true_sol, params] = step1(Target, Model, params, 3, Dataset);
            
            % Perform estimation and update parameter values
            [est_sol,  params] = step2(Target, Model, params, 3);
            
            % Run simulation using updated parameters
            [pred_sol, params] = step3(Target, Model, params, 3, est_sol);
    
            % Compute RMSE
            % Compute difference between data and simulation
            params = compute_RMSE(Target, params, baseline_sol, true_sol);
            % Unpack RMSE Values
            RMSE_values(i, j) = params.RMSE_mV;
            disp(params.RMSE_mV)

            % RMSE_values(i, j) = return_RMSE('Compare', varied_params, baseline_sol, pred_sol);
        end
    end
    
    % Create a table with delta values and corresponding RMSE
    delta_str = arrayfun(@(d) sprintf('Delta=%.2f', d), delta, 'UniformOutput', false);
    sensitivity_table = array2table(RMSE_values, 'RowNames', param_names_to_vary, 'VariableNames', delta_str);
    
    plot(delta, )
    % solution structure
    % compare true_sol with real data




    % Display results
    disp(sensitivity_table);
end
