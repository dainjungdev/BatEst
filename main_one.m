function [out, params] = main_one(Dataset,out,input_params,cycle_step,outputPath)
% This is the main script for a single simulation or optimisation step.
% The inputs and outputs are optional.

% close all;
reset_path;
startTime = datetime('now');
fprintf('\nCode started at %s\n', startTime);
% fprintf('Cycle: %d, Step: %d\n', cycle_step(1), cycle_step(2));

% Initialise optional variables
if ~exist('Dataset','var'), Dataset = import_parquet('Cell3_RPT0.parquet'); end
if ~exist('out','var'), out = []; end
if ~exist('input_params','var'), input_params = []; end
if ~exist('cycle_step','var'), cycle_step = []; end
if ~exist('outputPath', 'var'), outputPath = './Output'; end

%% Setup
% The following settings must be defined.
% ModelName: choose from the available Models (OCV, RORC, EHMT, etc.)
% Target: choose from Simulate, Plot, Compare or Parameter
% Estimator: choose from the available Methods (Fmincon, PEM)

% Settings
ModelName = 'OCV';
Target = 'Parameter';
Estimator = 'PEM';
DataType = 'Pseudo-OCV charge';
cycle_step =[0;10];

%% Start
fprintf('\nComputation started at %s\n', datetime("now"));

% Add relevant paths
reset_path;
addpath(genpath(strcat('./BatEst/Code/Models/',ModelName)));
addpath(genpath(strcat('./BatEst/Code/Methods/',Estimator)));

% Define dimensionless model
[Model, params] = step0(ModelName,1,input_params);
Model.Noise = false; % true or false
params.cycle_step = cycle_step;
params.DataType = DataType;

% Load or generate data
[true_sol, params] = step1(Target,Model,params,1,Dataset);

% Perform estimation and update parameter values
[est_sol,  params] = step2(Target,Model,params,1);

% Run simulation using updated parameters
[pred_sol, params] = step3(Target,Model,params,1,est_sol);

% Compare prediction and data
params = step4(Target,params,true_sol,pred_sol);

%% Save
% Only need to save the updated parameters structure as plots
% can be re-generated using Simulate or Compare.

% Convert and save the parameters in a table
out = tabulate_output(params, out);
% summary_table = summarise(out);

endTime = datetime('now');
duration = endTime - startTime;
fprintf('\nComputation ended at %s\n', endTime);
fprintf('Total duration: %s\n', duration);

% Save output and current figure
% outputPath = 'Project/06_cycling/out/cycle_OCV_test/Cell3_1_48_prev_estimate';
% mkdir(outputPath);
fileName = datestr(datetime('now'), 'yyyy-mm-dd_HH-MM-SS');
save_output(out, [outputPath '/' fileName], true);

% Save Summary
% writetable(summary_table, [outputPath '/' fileName '.csv']);
% save_plot(gcf,[outputPath '/' fileName]);

% reset_path;

end

