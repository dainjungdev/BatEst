function restoreSelectDataSubset()
    % Restore the original select_data_subset.m from backup
    backupFilename = 'BatEst/Code/Models/select_data_subset_backup.m';
    originalCode = fileread(backupFilename);
    fileID = fopen('BatEst/Code/Common/select_data_subset.m', 'w');
    fprintf(fileID, '%s', originalCode);
    fclose(fileID);
end
