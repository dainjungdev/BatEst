function [out, params] = main_multi(Dataset, out, params, rep_num, outputPath)

startTime = setupEnvironment();

% Initialise optional variables
if ~exist('Dataset','var'), Dataset = import_parquet('./BatEst/Data/Examples/Raj2020_Tests.parquet'); end
if ~exist('out','var'), out = []; end
if ~exist('params','var'), params = []; end
if ~exist('rep_num','var'), rep_num = 1:3; end
if ~exist('outputPath', 'var'), outputPath = './BatEst/Data/out/'; end

Dataset = [];

% Set the section number(s) or number of repetitions
for j = rep_num
    % Define Settings
    if j == 1
        ModelName = 'OCV';
    elseif j >= 2
        ModelName = 'EHM';
    end
    Target = 'Simulate';
    Estimator = 'PEM';

    addModelAndMethodPaths(ModelName, Estimator);
    
    [Model, params] = step0(ModelName, j, params);
    [true_sol, params] = step1(Target, Model, params, j, Dataset);
    [est_sol, params] = step2(Target, Model, params, j);
    [pred_sol, params] = step3(Target, Model, params, j, est_sol);
    params = step4(Target, params, true_sol, pred_sol);

    % Optionally save results
    out = tabulate_output(params, out);
end

finalizeComputation(startTime, outputPath, ModelName, out);
end

