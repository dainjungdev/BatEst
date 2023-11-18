function RMSE_values = sensitivity_real_data(input_params)
% Sensitivity analysis by comparing with real data

close all;
reset_path;

startTime = datetime('now');
fprintf('\nComputation started at %s\n', startTime);

% Load the dataset
Dataset = parquetread('Data/Examples/Raj2020_Tests.parquet');

% Define baseline parameters from input
baseline_params = input_params;
baseline_params.verbose = false;
baseline_params.plot_results = false;

% Model configuration
ModelName = 'EHM';
Estimator = 'PEM';
Target = 'Compare';
addpath(genpath(strcat('./Code/Models/', ModelName)));
addpath(genpath(strcat('./Code/Methods/', Estimator)));

% Run baseline simulation
[Model, baseline_params] = step0(ModelName, 3, baseline_params);
baseline_sol = step1(Target, Model, baseline_params, 3, Dataset);

% Define parameters to vary
varying_params = {'In', 'Ip', 'Q', 'Rf', 'tau', 'nu', 'miu', 'alph', 'b'};
num_parameters = length(varying_params);

% Define the range of delta
min_value = -0.15;
max_value = 0.15;
num_delta = 50;
delta = linspace(min_value, max_value, num_delta);

% Initialise RMSE results array
RMSE_values = zeros(num_parameters, length(delta));

% Loop over each parameter and delta
for i = 1:num_parameters
    param1 = varying_params{i};
    for j = 1:num_delta
        % Modify the parameter value
        varied_params = baseline_params;
        varied_params.verbose = false;
        varied_params.plot_results = false;

        varied_params.(param1) = baseline_params.(param1) * (1 + delta(j));

        % Run model simulation
        [Model, varied_params] = step0(ModelName, 3, varied_params);
        [true_sol, varied_params] = step1('Simulate', Model, varied_params, 3, Dataset);
        varied_params = step4('Compare', varied_params, baseline_sol, true_sol);

        % Store RMSE
        RMSE_values(i, j) = varied_params.RMSE_mV;
    end
end

% Plotting the results
figure; hold on;
for i = 1:num_parameters
    plot(delta, RMSE_values(i, :), 'DisplayName', varying_params{i});
end

xlabel('Parameter Variation Delta');
ylabel('RMSE (mV)');
title('Sensitivity Analysis for EHM Parameters');
legend('show');
hold off;

endTime = datetime('now');
duration = endTime - startTime;
fprintf('\nComputation ended at %s\n', endTime);
fprintf('Total duration: %s\n', duration);

% Adjust layout
set(gcf, 'Position', get(0, 'Screensize'));

% Save plot
save_plot(gcf,['1_sensitivity_analysis/out_' ModelName '_3_1'],true);
end
