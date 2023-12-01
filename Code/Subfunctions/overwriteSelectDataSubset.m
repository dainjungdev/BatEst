function overwriteSelectDataSubset(newCode)
    % Overwrite select_data_subset.m with new code

    fileID = fopen('BatEst/Code/Common/select_data_subset.m', 'w');
    fprintf(fileID, '%s', newCode);
    fclose(fileID);
end
