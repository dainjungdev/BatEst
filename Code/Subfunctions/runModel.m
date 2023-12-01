function params = runModel(ModelName, Target, Estimator, j, params, Dataset)
    addModelAndMethodPaths(ModelName, Estimator);
    
    [Model, params] = step0(ModelName, j, params);
    [true_sol, params] = step1(Target, Model, params, j, Dataset);
    [est_sol, params] = step2(Target, Model, params, j);
    [pred_sol, params] = step3(Target, Model, params, j, est_sol);
    params = step4(Target, params, true_sol, pred_sol);

    % Optionally save results
    out = tabulate_output(params, out);
end