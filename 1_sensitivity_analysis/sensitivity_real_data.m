function RMSE_values = sensitivity_real_data(input_params)
% Sensitivity analysis by comparing with real data

close all;
reset_path;

% Load the dataset
Dataset = parquetread('Data/Examples/Raj2020_Tests.parquet');

% Define baseline parameters from input
baseline_params = input_params;
baseline_params.verbose = false;

% Model configuration
ModelName = 'EHM';
Estimator = 'PEM';
Target = 'Compare';
addpath(genpath(strcat('./Code/Models/', ModelName)));
addpath(genpath(strcat('./Code/Methods/', Estimator)));

% Run baseline simulation
[Model, baseline_params] = step0(ModelName, 3, baseline_params);
baseline_sol = step1(Target, Model, baseline_params, 3, Dataset);

% Parameters to vary and their range
param_names_to_vary = {'In', 'Ip', 'Q', 'Rf', 'tau', 'nu', 'miu', 'alph', 'b'};
min_value = -0.15;
max_value = 0.15;
num_values = 50;
delta = linspace(min_value, max_value, num_values);

% Initialize RMSE results array
RMSE_values = zeros(length(param_names_to_vary), length(delta));

% Loop over each parameter and delta
for i = 1:length(param_names_to_vary)
    for j = 1:length(delta)
        % Modify the parameter value
        varied_params = baseline_params;
        varied_params.(param_names_to_vary{i}) = baseline_params.(param_names_to_vary{i}) * (1 + delta(j));
        varied_params.verbose = false;
        varied_params.plot_results = false;

        % Run simulation with varied parameters
        [Model, varied_params] = step0(ModelName, 3, varied_params);
        [true_sol, varied_params] = step1('Simulate', Model, varied_params, 3, Dataset);

        % Compare with baseline to compute RMSE
        varied_params = step4('Compare', varied_params, baseline_sol, true_sol);
        RMSE_values(i, j) = varied_params.RMSE_mV;
    end
end

% Plotting the results
figure; hold on;
for i = 1:length(param_names_to_vary)
    plot(delta, RMSE_values(i, :), 'DisplayName', param_names_to_vary{i});
end

xlabel('Parameter Variation Delta');
ylabel('RMSE (mV)');
title('Sensitivity Analysis for EHM Parameters');
legend('show');
hold off;
end
