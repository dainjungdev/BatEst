function [out, params] = main_multi_cycle(Dataset,out,input_params)
% Function similar to main_multi, for cycling data

close all;
reset_path;

startTime = datetime('now');
fprintf('\nmain_multi started at %s\n', startTime);

% Initialise optional variables
% To create a Dataset use: Dataset = importParquet('XXX.csv');
% Dataset = parquetread('Code/Common/Import/ExampleDataset.parquet');
if ~exist('Dataset','var'), Dataset = []; end

% To append the output from a previous estimation, use the input out.
if ~exist('out','var'), out = []; end

% To pass the parameters from a previous estimation, use the input params.
% To generate a parameters structure, use: params = load_output(out);
% input_params = load_output(out);
if ~exist('input_params','var'), input_params = []; end

cycles_data = import_parquet('Data/Examples/Raj2020_Cycling.parquet');

%% 2. Unique cycle indices
params = input_params;
cycle_data = cycles_data(cycles_data.Cycle_Index == 1, :);

% Settings
ModelName = 'EHM';
Target = 'Parameter';
Estimator = 'PEM';
Dataset = cycle_data;

%% Start
fprintf('\nComputation started at %s\n', datetime("now"));

% Add relevant paths
reset_path;
addpath(genpath(strcat('./Code/Models/',ModelName)));
addpath(genpath(strcat('./Code/Methods/',Estimator)));

% Define dimensionless model
[Model, params] = step0(ModelName,0,input_params);
Model.Noise = false; % true or false

% Load or generate data
[true_sol, params] = step1(Target,Model,params,6,Dataset);

% Perform estimation and update parameter values
[est_sol,  params] = step2(Target,Model,params,0);

% Run simulation using updated parameters
[pred_sol, params] = step3(Target,Model,params,0,est_sol);

% Compare prediction and data
params = step4(Target,params,true_sol,pred_sol);


%% Save
% Only need to save the updated parameters structure as plots
% can be re-generated using Simulate or Compare.

% Convert and save the parameters in a table
out = tabulate_output(params,out);

% Save output and current figure
save_output(out,['Data/out_' ModelName]);
% save_plot(gcf,['Data/plot_' ModelName]);


end