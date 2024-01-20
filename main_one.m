function [out, params] = main_one(Dataset,out,input_params,outputPath)
% This is the main script for a single simulation or optimisation step.
% The inputs and outputs are optional.

startTime = setupEnvironment();

% Initialise optional variables
if ~exist('Dataset','var'), Dataset = []; end
if ~exist('out','var'), out = []; end
if ~exist('input_params','var'), input_params = []; end
if ~exist('outputPath', 'var'), outputPath = './BatEst/Data/out'; end


%% Setup
% The following settings must be defined.
% ModelName: choose from the available Models (OCV, RORC, EHMT, etc.)
% Target: choose from Simulate, Plot, Compare or Parameter
% Estimator: choose from the available Methods (Fmincon, PEM)

% Settings
ModelName = 'OCV_MSMR';
Target = 'Simulate';
Estimator = 'PEM';


%% Start
fprintf('\nComputation started at %s\n', datetime("now"));

% Add relevant paths
reset_path;
addpath(genpath(strcat('./BatEst/Code/Models/',ModelName)));
addpath(genpath(strcat('./BatEst/Code/Methods/',Estimator)));

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

out = tabulate_output(params, out);

finalizeComputation(startTime, outputPath, ModelName, out);

end

