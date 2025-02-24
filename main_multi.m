function [out, params] = main_multi(Dataset, out, params, cycle, outputPath)

close all;
reset_path;
startTime = datetime('now');
fprintf('\nCode started at %s\n', startTime);

% Initialize optional variables
if ~exist('Dataset','var'), Dataset = import_parquet('Cell4_RPT0.parquet'); end
if ~exist('out','var'), out = []; end
if ~exist('params','var'), params = []; end
if ~exist('cycle','var'), cycle = 0*48; end
if ~exist('outputPath', 'var'), outputPath = './Output'; end


fprintf('\n-------\nCycle: %d\n', cycle);
% params = load("Cell20_Cycle96_Step3.mat").params;

rep_num = 1:3;
% Set the section number(s) or number of repetitions
for j = rep_num
    % Define Settings
    if j == 1
        ModelName = 'OCV';
        cycle_step = [cycle;10];
        DataType = 'Pseudo-OCV charge';
    elseif j == 2
        ModelName = 'EHM';
        cycle_step = [cycle;5];
        DataType = 'Relaxation';
    elseif j == 3
        ModelName = 'EHM';
        cycle_step = [cycle;6];
        DataType = 'CCCV charge';
    end

    Target = 'Parameter';
    Estimator = 'PEM';

    addpath(genpath(strcat('./BatEst/Code/Models/', ModelName)));
    addpath(genpath(strcat('./BatEst/Code/Methods/', Estimator)));

    % Define dimensionless model
    [Model, params] = step0(ModelName,j,params);
    Model.Noise = false; % true or false
    params.DataType = DataType;
    params.cycle_step = cycle_step;
    
    % Load or generate data
    [true_sol, params] = step1(Target,Model,params,j,Dataset);

    % Perform estimation and update parameter values
    [est_sol,  params] = step2(Target,Model,params,j);
    
    % Run simulation using updated parameters
    [pred_sol, params] = step3(Target,Model,params,j,est_sol);
    
    % Compare prediction and data
    params = step4(Target,params,true_sol,pred_sol);
    
    % Optionally save results
    out = tabulate_output(params, out);
    summary_table = summarise(out);
end

%% Save
% Only need to save the updated parameters structure as plots
% can be re-generated using Simulate or Compare.

endTime = datetime('now');
duration = endTime - startTime;
fprintf('\nComputation ended at %s\n', endTime);
fprintf('Total duration: %s\n', duration);

% Save output and current figure
% outputPath = 'Project/06_cycling/out/cycle_OCV_test/Cell3_1_48_prev_estimate';
% mkdir(outputPath);
fileName = datestr(datetime('now'), 'yyyy-mm-dd_HH-MM-SS');
% save(fileName, "params");
% save_output(out, [outputPath '/' fileName], true);
save_output(out, [outputPath '/' fileName '_' num2str(cycle)], true);

% Save Summary
% writetable(summary_table, [outputPath '/' fileName '.csv']);
writetable(summary_table, [outputPath '/' fileName '_' num2str(cycle) '.csv']);
% save_plot(gcf,[outputPath '/' fileName]);

% reset_path;
end

