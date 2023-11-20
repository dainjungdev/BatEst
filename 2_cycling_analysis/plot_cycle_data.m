function plot_cycle_data(file_path, is_subplot)
% This function plots current and voltage over time for each cycle
% is_subplot: if true, plot as subplots, else separate plots

    % Load data
    % data = parquetread(file_path);
    data = parquetread('Data/Examples/Raj2020_Cycling.parquet');
    
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
