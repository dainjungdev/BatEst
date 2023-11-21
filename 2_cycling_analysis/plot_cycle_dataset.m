function plot_cycle_dataset(file_path, is_subplot)
% This function plots current and voltage over time for each cycle
% is_subplot: if true, plot as subplots, else separate plots

close all;
reset_path;

startTime = datetime('now');
fprintf('\nComputation started at %s\n', startTime);

% Load dataset
% data = parquetread(file_path);
% When 'file_path' does
if ~exist('file_path','var'), 
    inputStruct = load('2_cycling_analysis/TPG2_data/TPG2.2-Cell3.mat');
    fieldName = fieldnames(inputStruct);
    data = inputStruct.(fieldName{1});
else
    folderPath = file_path;
end

if ~exist('is_subplot','var'), 
    is_subplot = true; 
end

% List all .mat files in the folder
matFiles = dir(fullfile(folderPath, '*.mat'));

for j = 1:length(matFiles)
    % Load the table from the .mat file
    filePath = fullfile(folderPath, matFiles(j).name);
    loadedData = load(filePath);
    fieldName = fieldnames(loadedData);
    data = loadedData.(fieldName{1});

    % Unique cycle indices
    cycles = unique(data.Cycle_Index);

    % Preparing color map for different cycles
    colors = turbo(length(cycles));  % Or choose another colormap

    %% Create a single figure with subplots
    if is_subplot
        figure('Name', 'Data Analysis - Subplots');
    
        % Subplot for Current
        subplot(2, 1, 1);
        hold on;
        for i = 1:length(cycles)
            cycle_data = data(data.Cycle_Index == cycles(i), :);
            plot(cycle_data.Test_Time_s, cycle_data.Current_A, '.', 'Color', colors(i, :));
        end
        xlabel('Time (s)');
        ylabel('Current (A)');
        title('Current over Time for Each Cycle');
        legend(string(cycles));
        hold off;
    
        % Subplot for Voltage
        subplot(2, 1, 2);
        hold on;
        for i = 1:length(cycles)
            cycle_data = data(data.Cycle_Index == cycles(i), :);
            plot(cycle_data.Test_Time_s, cycle_data.Voltage_V, '.', 'Color', colors(i, :));
        end
        xlabel('Time (s)');
        ylabel('Voltage (V)');
        title('Voltage over Time for Each Cycle');
        legend(string(cycles));
        hold off;
    
        % Adjust layout
        set(gcf, 'Position', get(0, 'Screensize'));
        
        % Save plot
        save_plot(gcf,['2_cycling_analysis/out/plot_cycle_dataset/' fieldName{1}],true);
        endTime = datetime('now');
        duration = endTime - startTime;
        fprintf('\nComputation ended at %s\n', endTime);
        fprintf('Total duration: %s\n', duration);
    
    %% Creating separate figures for each plot
    else
        % Figure for Current plot
        figure('Name', 'Current over Time for Each Cycle');
        hold on;
        for i = 1:length(cycles)
            cycle_data = data(data.Cycle_Index == cycles(i), :);
            plot(cycle_data.Test_Time_s, cycle_data.Current_A, '.', 'Color', colors(i, :));
        end
        xlabel('Time (s)');
        ylabel('Current (A)');
        title('Current over Time for Each Cycle');
        legend(string(cycles));
        hold off;
    
        % Figure for Voltage plot
        figure('Name', 'Voltage over Time for Each Cycle');
        hold on;
        for i = 1:length(cycles)
            cycle_data = data(data.Cycle_Index == cycles(i), :);
            plot(cycle_data.Test_Time_s, cycle_data.Voltage_V, '.', 'Color', colors(i, :));
        end
        xlabel('Time (s)');
        ylabel('Voltage (V)');
        title('Voltage over Time for Each Cycle');
        legend(string(cycles));
        hold off;
    end    
end
end
