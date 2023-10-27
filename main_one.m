function [out, params] = main_one(Dataset, out, input_params)
% main_one: main script for a single simulation or optimisation step
% The inputs and outputs are optional.

%% 0. Initial Step
close all;
reset_path;


%% 1. Initialise optional variables
% Create dataset(parquet) from csv: Dataset = importParquet('XXX.csv')
if ~exist('Dataset', 'var')
    Dataset = [];
end

% Append output from previous estimation: out
if ~exist('out', 'var')
    out = [];
end

% Generate parameters from previous estimation: input_params = load_output(out)
if ~exist('input_params', 'var')
    input_params = [];
end


%% 2. Setup
% The following settings must be defined.
% ModelName: choose from the available Models (OCV, RORC, EHMT, etc.)
% Target: choose from Simulate, Plot, Compare or Parameter
% Estimator: choose from the available Methods (Fmincon, PEM)

ModelName = 'EHMT';
Target = 'Simulate';
Estimator = 'PEM';

% Add relevant paths
reset_path;
addpath(genpath(strcat('./Code/Models/', ModelName)));
addpath(genpath(strcat('./Code/Methods/', Estimator)));


%% 3. Computation
fprintf('\nComputation started at %s\n', datetime("now"));

% Define dimensionless model
[Model, params] = step0(ModelName, 0, input_params);
Model.Noise = false; % true or false

% Load or generate data
[true_sol, params] = step1(Target, Model, params, 6, Dataset);

% Perform estimation and update parameter values
[est_sol, params] = step2(Target, Model, params, 0);

% Run simulation using updated parameters
[pred_sol, params] = step3(Target, Model, params, 0, est_sol);

% Compute difference between data and simulation
params = compute_RMSE(Target, params, true_sol, pred_sol);


%% 4. Save
% Only need to save the updated parameters structure
% Plots can be re-generated using Simulate or Compare.

% Convert and save the parameters in a table
out = tabulate_output(params, out);

% Save output and current figure
save_output(out,['Data/out_' ModelName]);
save_plot(gcf,['Data/plot_' ModelName]);
end

