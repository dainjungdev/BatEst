function [ModelName, Target, Estimator] = defineSettings(j)
    if j == 1
        ModelName = 'OCV';
    elseif j >= 2
        ModelName = 'EHM';
    end
    Target = 'Parameter';
    Estimator = 'PEM';
    % Update Dataset if needed
    % Dataset remains unchanged if not specified here
end