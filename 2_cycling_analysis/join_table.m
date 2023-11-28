function joined_table = join_table(folderPath)
    % List all .parquet files in the folder
    parquetFiles = dir(fullfile(folderPath, '*.parquet'));

    filePath = fullfile(folderPath, parquetFiles(1).name);
    joined_table = parquetread(filePath);

    for i = 2:length(parquetFiles)
        % Load the .mat file
        filePath = fullfile(folderPath, parquetFiles(i).name);
        loadedData = parquetread(filePath);
        joined_table(i,:) = loadedData(1,:);
    end