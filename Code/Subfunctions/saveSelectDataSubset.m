function saveSelectDataSubset(outputPath, selectDataSubsetPath)
    selectDataSubsetContent = fileread(selectDataSubsetPath);
    savedFilename = fullfile(outputPath, 'data_selection.txt');
    fileID = fopen(savedFilename, 'w');
    fprintf(fileID, '%s', selectDataSubsetContent);
    fclose(fileID);
end