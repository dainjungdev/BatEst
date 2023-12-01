function [cycle_step, DataType] = select_data_subset(j)
% Select subset of data based on index j
switch j
    case 1
        cycle_step = [0;10];
        DataType = 'Pseudo-OCV charge';
    case 2
        cycle_step = [0;5];
        DataType = 'Relaxation';
    case 3
        cycle_step = [0;6];
        DataType = 'CCCV charge';
    otherwise
        % Create empty properties
        cycle_step = [];
        DataType = '';
end
end