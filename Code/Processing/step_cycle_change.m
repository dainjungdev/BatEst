function changeIndices = step_cycle_change(myTable)
    % Get the arrays of Step_Index and Cycle_Index values
    stepIndexValues = myTable.Step_Index;
    cycleIndexValues = myTable.Cycle_Index;

    % Find indices where Step_Index changes
    stepChangeRows = find([true; diff(stepIndexValues) ~= 0]);

    % Get the corresponding Step_Index values
    stepIndexChanges = stepIndexValues(stepChangeRows);

    % Find indices where Cycle_Index changes
    cycleChangeRows = find([true; diff(cycleIndexValues) ~= 0]);

    % Get the corresponding Cycle_Index values
    cycleIndexChanges = cycleIndexValues(cycleChangeRows);

    % Combine row indices and Step_Index values
    changeIndices = [stepChangeRows, stepIndexChanges, cycleIndexValues(stepChangeRows)];

    % If there's a change in Cycle_Index at the same row, add it
    for i = 1:length(cycleChangeRows)
        rowIndex = cycleChangeRows(i);
        cycleChange = cycleIndexChanges(i);

        % Check if the change in Cycle_Index corresponds to an existing row
        if ismember(rowIndex, changeIndices(:, 1))
            changeIndices(changeIndices(:, 1) == rowIndex, 3) = cycleChange;
        else
            changeIndices = [changeIndices; rowIndex, stepIndexValues(rowIndex), cycleChange];
        end
    end

    % Sort the changeIndices by row number
    changeIndices = sortrows(changeIndices, 1);
end
