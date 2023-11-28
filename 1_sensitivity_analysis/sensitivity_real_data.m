function [MAE_table, RMSE_table, MaxError_table] = sensitivity_real_data(input_params)
% Sensitivity analysis by comparing with real data
% L1 Norm: Mean Absolute Error
% L2 Norm: Root Mean Square Error
% L-infinity Norm : Maximum Absolute Error

close all;
reset_path;

startTime = datetime('now');
fprintf('\nComputation started at %s\n', startTime);

% Load dataset
Dataset = parquetread('Data/Examples/Raj2020_Tests.parquet');

% Define baseline parameters
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
varying_params = {'Q', 'tau', 'b', 'Ip', 'In', 'nu', 'miu', 'Rf'};
num_parameters = length(varying_params);

% Define the range of delta
min_value = -0.1;
max_value = 0.1;
num_delta = 11; % odd number
delta = linspace(min_value, max_value, num_delta);


%% Compute RMSE for each parameter
% Initialise results array
MAE_values = zeros(num_parameters, length(delta)); % L1 norm
RMSE_values = zeros(num_parameters, length(delta)); % L2 norm
MaxError_values = zeros(num_parameters, length(delta)); % L-infinity norm

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


        % Store MAE
        MAE_values(i, j) = varied_params.MAE_mV;

        % Store RMSE
        RMSE_values(i, j) = varied_params.RMSE_mV;

        % Store Max Absolute Error
        MaxError_values(i, j) = varied_params.MaxError_mV;
    end
end


%% Compute Standard Deviation of RMSE for each parameter
std_RMSE = std(RMSE_values, 0, 2);

% Convert RMSE_values to a table
RMSE_table = array2table(RMSE_values, 'VariableNames', strcat('Delta=', string(delta)));
RMSE_table.Properties.RowNames = varying_params;
RMSE_table.std_RMSE = std_RMSE;
save_output(RMSE_table,['1_sensitivity_analysis/new' ...
    num2str(min_value) '_' num2str(max_value) '_' num2str(num_delta)], true)

if num_delta < 10
    disp(RMSE_table);
else
    std_RMSE_struct = struct();
    for i = 1:num_parameters
        std_RMSE_struct.(varying_params{i}) = std_RMSE(i);
    end
    
    disp('Standard Deviation of RMSE for each parameter:');
    disp(std_RMSE_struct);
end

% Convert MAE_values and MaxError_values to tables
MAE_table = array2table(MAE_values, 'VariableNames', strcat('Delta=', string(delta)));
MAE_table.Properties.RowNames = varying_params;

MaxError_table = array2table(MaxError_values, 'VariableNames', strcat('Delta=', string(delta)));
MaxError_table.Properties.RowNames = varying_params;


%% Plotting the results : RMSE
figure; hold on;
for i = 1:num_parameters
    plot(delta, RMSE_values(i, :), 'DisplayName', varying_params{i});
end

% Calculate interquartile range (IQR) and determine y-axis limits
all_RMSE_values = RMSE_values(:); % Convert 2D array to a vector
Q1 = quantile(all_RMSE_values, 0.25);
Q3 = quantile(all_RMSE_values, 0.75);
IQR = Q3 - Q1;

% Set y-axis limits to exclude outliers
lowerLimit = max(min(all_RMSE_values), Q1 - 2 * IQR);
upperLimit = min(max(all_RMSE_values), Q3 + 2 * IQR);
ylim([lowerLimit upperLimit]);

xlabel('Parameter Variation Delta');
ylabel('RMSE (mV)');
mainTitle = 'Sensitivity Analysis for EHM Parameters';
deltaInfo = sprintf('Delta range: [%.2f, %.2f], Number of Deltas: %d', min_value, max_value, num_delta);
title({mainTitle, deltaInfo});  % Use a cell array for multi-line title

legend('show');
hold off;

% Adjust layout
set(gcf, 'Position', get(0, 'Screensize'));

% Save plot
save_plot(gcf,['1_sensitivity_analysis/out/sensitivity_real_' num2str(min_value) '_' num2str(max_value) '_' num2str(num_delta)],true);
endTime = datetime('now');
duration = endTime - startTime;
fprintf('\nComputation ended at %s\n', endTime);
fprintf('Total duration: %s\n', duration);

%% Additional plotting for MAE and Max Absolute Error
figure; 
subplot(2,1,1); hold on;
for i = 1:num_parameters
    plot(delta, MAE_values(i, :), 'DisplayName', varying_params{i});
end
xlabel('Parameter Variation Delta');
ylabel('MAE (mV)');
title('MAE Sensitivity Analysis for EHM Parameters');
legend('show');
hold off;

subplot(2,1,2); hold on;
for i = 1:num_parameters
    plot(delta, MaxError_values(i, :), 'DisplayName', varying_params{i});
end
xlabel('Parameter Variation Delta');
ylabel('Max Absolute Error (mV)');
title('Max Absolute Error Sensitivity Analysis for EHM Parameters');
legend('show');
hold off;

% Adjust layout
set(gcf, 'Position', get(0, 'Screensize'));

% Save additional metrics tables
save_output(MAE_table,['1_sensitivity_analysis/out_sensitivity_MAE_table' ...
    num2str(min_value) '_' num2str(max_value) '_' num2str(num_delta)], true);
save_output(MaxError_table,['1_sensitivity_analysis/out_sensitivity_MaxError_table' ...
    num2str(min_value) '_' num2str(max_value) '_' num2str(num_delta)], true);

end