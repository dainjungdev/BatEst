function addModelAndMethodPaths(ModelName, Estimator)
    addpath(genpath(strcat('./BatEst/Code/Models/', ModelName)));
    addpath(genpath(strcat('./BatEst/Code/Methods/', Estimator)));
end