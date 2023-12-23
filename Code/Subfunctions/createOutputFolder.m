function outputPath = createOutputFolder(basePath)
    folderName = datestr(datetime('now'), 'yyyy-mm-dd_HH-MM-SS');
    outputPath = [basePath folderName];
    mkdir(outputPath);
end