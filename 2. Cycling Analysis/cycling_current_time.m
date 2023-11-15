data = parquetread('Data/Examples/Raj2020_Cycling.parquet');

% Unique cycle indices
cycles = unique(data.Cycle_Index);

% Preparing color map for different cycles
colors = turbo(length(cycles));  % Or choose another colormap

% Creating figure for Current plot
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

% Creating figure for Voltage plot
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
