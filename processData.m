function processData(folderPath)
% Changes all files in folder(Raj. Dataset) to suitable format for BatEst
    % List all .mat files in the folder
    matFiles = dir(fullfile(folderPath, '*.mat'));

    % Create 'out' directory if it doesn't exist
    outFolder = fullfile(folderPath, 'out');
    if ~exist(outFolder, 'dir')
        mkdir(outFolder);
    end

    % Process each .mat file
    for i = 1:length(matFiles)
        % Load the table from the .mat file
        filePath = fullfile(folderPath, matFiles(i).name);
        loadedData = load(filePath);
        tableName = fieldnames(loadedData);
        originalTable = loadedData.(tableName{1});

        % Process the table
        processedTable = processTable(originalTable); 

        % Save the processed table in the 'out' folder with the same filename
        saveFilePath = fullfile(outFolder, matFiles(i).name);

        % Create a structure with a dynamic field name
        S.(tableName{1}) = processedTable;

        % Save the structure in the .mat file
        save(saveFilePath, '-struct', 'S');

        % Clear the temporary structure for the next iteration
        clear S;
    end
end
