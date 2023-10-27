function [out, params] = main_multi(Dataset, out, input_params)
% main_multi: main script to run multiple iterations
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


%% 2. Iterations
% Load the index of measurement data
index_name = 'Data/Examples/Test_Index.parquet';
index = parquetread(index_name);

% Set the cell numbers to process
cell_nums = [3];

% Iterate through each cell number
for cell_num = cell_nums
    
    % Initialise parameters for the current cell
    params = input_params;
    
    % Find relevant files based on cell number
    file_indices = find(index.Cell_Number == cell_num & index.Performance_Test);
    filenames = index.File_Name(file_indices);
    subfolders = index.Folder_Name(file_indices);
    
    % Loop through each file for the current cell
    for file_idx = 1:length(filenames)
        
        % Initialise parameters for the current file
        params = input_params;
        
        % Set section numbers or number of repetitions
        rep_nums = 1:3;
        
        % Iterate through each repetition or section
        for rep_num = rep_nums

            %% 2-1. Setup
            % The following settings must be defined.
            % ModelName: choose from the available Models (OCV, RORC, EHMT, etc.)
            % Target: choose from Simulate, Plot, Compare or Parameter
            % Estimator: choose from the available Methods (Fmincon, PEM)

            % Initialise settings based on repetition number
            if rep_num == 1
                ModelName = 'OCV';  % For pseudo-OCV measurements
                Target = 'Compare';
                Estimator = 'PEM';
                Dataset = import_parquet([subfolders{file_idx} '/' filenames{file_idx}]);
            elseif rep_num == 2
                ModelName = 'EHM';  % For relaxation data
                % ... (Other settings)
            elseif rep_num == 3
                ModelName = 'EHM';  % For CCCV charge data
                % ... (Other settings)
            end
            
            % Add relevant paths
            reset_path;
            addpath(genpath(strcat('./Code/Models/',ModelName)));
            addpath(genpath(strcat('./Code/Methods/',Estimator)));


            %% 2-2. Computation
            fprintf('\nRepetition Number = %d\nModelName = %s\n', rep_num, ModelName);
            fprintf('\nComputation started at %s\n', datetime("now"));
            
            % Define dimensionless model
            [Model, params] = step0(ModelName, rep_num, input_params);
            Model.Noise = false; % true or false
            
            % Load or generate data
            [true_sol, params] = step1(Target, Model, params, rep_num, Dataset);
            
            % Perform estimation and update parameter values
            [est_sol,  params] = step2(Target, Model, params, rep_num);
            
            % Run simulation using updated parameters
            [pred_sol, params] = step3(Target, Model, params, rep_num, est_sol);
            
            % Compute difference between data and simulation
            params = compute_RMSE(Target, params, true_sol, pred_sol);


            %% 2-3. Save
            % Only need to save the updated parameters structure
            % Plots can be re-generated using Simulate or Compare.
            
            % Convert and save the parameters in a table
            out = tabulate_output(params, out);
            
            % Save output and current figure (true = overwrite by default)
            save_output(out, ['Data/out_' ModelName '_' num2str(cell_num) '_' num2str(file_idx)], true);
            % save_plot(gcf, ['Data/plot_' ModelName '_' num2str(cell_num) '_' num2str(file_idx)], true);
        end
    end
end

end

