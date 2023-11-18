function RMSE_matrix = sensitivity_real_data_2d(input_params)
% Sensitivity analysis by comparing with real data, combinations of two elements

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
% varying_params = {'Q', 'tau', 'b', 'Ip', 'In', 'nu', 'miu', 'Rf'};
varying_params = {'Q', 'In', 'nu', 'miu', 'Rf'};

% Get the indices of all possible combinations of two different elements
comb_indices = nchoosek(1:length(varying_params), 2);
num_combinations = size(comb_indices, 1);

% Define the range of delta
min_value = -0.05;
max_value = 0.05;
num_delta = 21; % odd number
delta = linspace(min_value, max_value, num_delta);


%% Compute RMSE for each parameter
% Initialise RMSE matrix
RMSE_matrix = zeros(num_combinations, num_delta, num_delta);

% Loop over each combination of parameters
for k = 1:num_combinations
    param1 = varying_params{comb_indices(k, 1)};
    param2 = varying_params{comb_indices(k, 2)};

    for i = 1:num_delta
        for j = 1:num_delta
            % Modify the parameter value
            varied_params = baseline_params;
            varied_params.verbose = false;
            varied_params.plot_results = false;

            varied_params.(param1) = baseline_params.(param1) * (1 + delta(i));
            varied_params.(param2) = baseline_params.(param2) * (1 + delta(j));
            
            % Run model simulation
            [Model, varied_params] = step0(ModelName, 3, varied_params);
            [true_sol, varied_params] = step1('Simulate', Model, varied_params, 3, Dataset);
            varied_params = step4('Compare', varied_params, baseline_sol, true_sol);
        
            % Store RMSE
            RMSE_matrix(k, i, j) = varied_params.RMSE_mV;
        end
    end
end



%% Plotting the results
figure;

% Calculate the number of subplots needed
num_subplots_x = ceil(sqrt(num_combinations));
num_subplots_y = ceil(num_combinations/num_subplots_x);

for k = 1:num_combinations
    param1 = varying_params{comb_indices(k, 1)};
    param2 = varying_params{comb_indices(k, 2)};
    
    % Generate a meshgrid for the parameter variations
    [Param1Grid, Param2Grid] = meshgrid(delta, delta);
    
    % RMSE values for the current combination
    CurrentRMSE = squeeze(RMSE_matrix(k, :, :));
    
    % Find the minimum RMSE value and its corresponding indices
    [minValue, minIndex] = min(CurrentRMSE(:));
    [minRow, minCol] = ind2sub(size(CurrentRMSE), minIndex);
    minParam1 = Param1Grid(minCol, minRow);
    minParam2 = Param2Grid(minCol, minRow);

    % Create 3D surface plot with a grid for each combination
    subplot(num_subplots_x, num_subplots_y, k);
    surf(Param1Grid, Param2Grid, CurrentRMSE', 'EdgeColor', 'k');
    colorbar; 
    xlabel(param1);
    ylabel(param2);
    zlabel('RMSE');
    title(sprintf('RMSE (%s vs %s)', param1, param2));
    
    view([-30 30]);
    hold on;
    grid on;
    
    % Highlight minimum RMSE value
    scatter3(minParam1, minParam2, minValue, 'filled', 'o', 'MarkerEdgeColor', ...
        'k', 'MarkerFaceColor', 'r', 'SizeData', 100);
    deltaParam1 = delta(minRow);
    deltaParam2 = delta(minCol);
    subtitle(sprintf('Min: %.3f, %s: %.3f, %s: %.3f', minValue, param1, deltaParam1, param2, deltaParam2))
    hold off;
end

% Adjust layout
set(gcf, 'Position', get(0, 'Screensize'));

% Save plot
save_plot(gcf,['1_sensitivity_analysis/out/sensitivity_real_2d' num2str(min_value) '_' num2str(max_value) '_' num2str(num_delta)],true);

endTime = datetime('now');
duration = endTime - startTime;
fprintf('\nComputation ended at %s\n', endTime);
fprintf('Total duration: %s\n', duration);
end
