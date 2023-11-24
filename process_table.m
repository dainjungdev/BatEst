function processedTable = processTable(originalTable)
% Table in Raj. Dataset to suitable format for BatEst
    % Rename columns
    originalTable = renamevars(originalTable, ...
                               {'Cyc', 'Step', 'TestTime', 'Amps', 'Volts', 'Temp1'}, ...
                               {'Cycle_Index', 'Step_Index', 'Test_Time_s', 'Current_A', 'Voltage_V', 'Temperature_C'});

    % Select only the required columns
    requiredColumns = {'Cycle_Index', 'Step_Index', 'Test_Time_s', 'Current_A', 'Voltage_V', 'Temperature_C'};
    processedTable = originalTable(:, requiredColumns);
end
