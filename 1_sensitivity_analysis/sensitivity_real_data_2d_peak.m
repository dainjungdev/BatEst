function sensitivity_real_data_2d_peak(input_params)
% Sensitivity analysis with real data but at specific data points

close all;
reset_path;

startTime = datetime('now');
fprintf('\nComputation started at %s\n', startTime);

% Load dataset
Dataset = parquetread('Data/Examples/Raj2020_Tests.parquet');

% Define baseline parameters
baseline_params = input_params;
% baseline_params.verbose = false;
% baseline_params.plot_results = false;

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
varying_params = {'nu', 'miu'};

% Set comb_indices manually
comb_indices = [1, 2];
num_combinations = size(comb_indices, 1);

% Set delta manually
% Define the range of delta
delta1 = -0.02;
delta2 = 0.03;

%% Compute RMSE for each parameter
% Loop over each combination of parameters
param1 = varying_params{comb_indices(1)};
param2 = varying_params{comb_indices(2)};

% Modify the parameter value
varied_params = baseline_params;
% varied_params.verbose = false;
% varied_params.plot_results = false;

varied_params.(param1) = baseline_params.(param1) * (1 + delta1);
varied_params.(param2) = baseline_params.(param2) * (1 + delta2);

% Run model simulation
[Model, varied_params] = step0(ModelName, 3, varied_params);
[true_sol, varied_params] = step1('Simulate', Model, varied_params, 3, Dataset);
varied_params = step4('Compare', varied_params, baseline_sol, true_sol);

disp(varied_params.RMSE_mV)

endTime = datetime('now');
duration = endTime - startTime;
fprintf('\nComputation ended at %s\n', endTime);
fprintf('Total duration: %s\n', duration);
end
