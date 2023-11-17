function RMSE_matrix = sensitivity_real_data_2d(input_params)
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

% % Parameters to vary and their range
% param_names_to_vary = {'In', 'Ip', 'Q', 'Rf', 'tau', 'nu', 'miu', 'alph', 'b'};
% min_value = -0.15;
% max_value = 0.15;
% num_values = 50;
% delta = linspace(min_value, max_value, num_values);

% Define ranges for two parameters
param_names_to_vary = {'Q', 'nu', 'miu'};
min1 = -0.15; max1 = 0.15; num_values1 = 5;
min2 = -0.05; max2 = 0.05; num_values2 = 5;
delta1 = linspace(min1, max1, num_values1);
delta2 = linspace(min2, max2, num_values2);

% % Initialize RMSE results array
% RMSE_values = zeros(length(param_names_to_vary), length(delta));

% Initialize RMSE matrix
RMSE_matrix = zeros(length(param_names_to_vary), length(param_names_to_vary), length(delta1), length(delta2));

% % Loop over each parameter and delta
% for i = 1:length(param_names_to_vary)
%     for j = 1:length(delta)
%         % Modify the parameter value
%         varied_params = baseline_params;
%         varied_params.(param_names_to_vary{i}) = baseline_params.(param_names_to_vary{i}) * (1 + delta(j));
%         varied_params.verbose = false;
%         varied_params.plot_results = false;
% 
%         % Run simulation with varied parameters
%         [Model, varied_params] = step0(ModelName, 3, varied_params);
%         [true_sol, varied_params] = step1('Simulate', Model, varied_params, 3, Dataset);
% 
%         % Compare with baseline to compute RMSE
%         varied_params = step4('Compare', varied_params, baseline_sol, true_sol);
%         RMSE_values(i, j) = varied_params.RMSE_mV;
%     end
% end

% Loop over each combination of parameters
for ii = 1:length(param_names_to_vary)
    for jj = ii:length(param_names_to_vary)
        for i = 1:length(delta1)
            for j = 1:length(delta2)
                % Set the parameter values
                param1 = param_names_to_vary{ii};
                param2 = param_names_to_vary{jj};

                varied_params = baseline_params;
                varied_params.verbose = false;
                varied_params.plot_results = false;

                varied_params.(param1) = baseline_params.(param1) * (1 + delta1(i));
                varied_params.(param2) = baseline_params.(param2) * (1 + delta2(j));

                % Run the model simulation
                [Model, varied_params] = step0(ModelName, 3, varied_params);
                [true_sol, varied_params] = step1('Simulate', Model, varied_params, 3, Dataset);
                varied_params = step4('Compare', varied_params, baseline_sol, true_sol);

                % Store RMSE
                RMSE_matrix(ii, jj, i, j) = varied_params.RMSE_mV;
            end
        end

        % Generate a meshgrid for the parameter ranges
        [Param1Grid, Param2Grid] = meshgrid(delta1, delta2);
        
        % RMSE values for the current combination of parameters
        CurrentRMSE = squeeze(RMSE_matrix(ii, jj, :, :));
        
        % Create the 3D surface plot
        figure;
        surf(Param1Grid, Param2Grid, CurrentRMSE');
        shading interp; % Optional: to smooth the colors on the surface
        colorbar;
        xlabel(param_names_to_vary{ii});
        ylabel(param_names_to_vary{jj});
        zlabel('RMSE');
        title('3D Surface of RMSE for Parameter Sensitivity');
        
        % Adjusting the view to match the example
        view([-30 30]); % Adjust the view angle for better visibility
        grid on;

        % Optionally, you can set the grid style
        set(gca, 'GridLineStyle', ':', 'GridColor', 'k', 'GridAlpha', 0.6); % Dotted grid lines
    end
end

% % Plotting the results
% figure; hold on;
% for i = 1:length(param_names_to_vary)
%     plot(delta, RMSE_values(i, :), 'DisplayName', param_names_to_vary{i});
% end
% 
% xlabel('Parameter Variation Delta');
% ylabel('RMSE (mV)');
% title('Sensitivity Analysis for EHM Parameters');
% legend('show');
% hold off;


end
