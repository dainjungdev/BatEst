function plot_cycle_compare(file_path)
% This function plots voltage over time for each cycle
% is_subplot: if true, plot as subplots, else separate plots

close all;
reset_path;

startTime = datetime('now');
fprintf('\nComputation started at %s\n', startTime);

% Load dataset
% data = parquetread(file_path);
if ~exist('file_path','var'), 
    data = parquetread('Data/Examples/Raj2020_Cycling.parquet'); 
end

% Unique cycle indices
% cycles = unique(data.Cycle_Index);
cycles = [1;2;3];

% Preparing color map for different cycles
colors = ['r', 'g', 'b'];

% Figure for Voltage plot
figure('Name', 'Voltage over Time for Each Cycle');
hold on;
for i = 1:length(cycles)
    cycle_data = data(data.Cycle_Index == cycles(i), :);
    % Subtract the first Test_Time_s value from all Test_Time_s values in the cycle
    relative_time = cycle_data.Test_Time_s - cycle_data.Test_Time_s(1);
    plot(relative_time, cycle_data.Voltage_V, '.', 'Color', colors(i));
end
xlabel('Relative Time (s)');
ylabel('Voltage (V)');
title('Voltage over Relative Time for Each Cycle');
legend(string(cycles));

% Adjust layout
set(gcf, 'Position', get(0, 'Screensize'));

% Save plot
save_plot(gcf,'2_cycling_analysis/out/plot_cycle_compare',true);

endTime = datetime('now');
duration = endTime - startTime;
fprintf('\nComputation ended at %s\n', endTime);
fprintf('Total duration: %s\n', duration);

end
