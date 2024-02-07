function [SimOpts, EstOpts] = DefaultPEMOpts
% Set the default simulation options for PEM.

SimOpts = struct('Solver','Auto', 'RelTol',1e-3, 'AbsTol',1e-5, ...
                 'MinStep','Auto', 'MaxStep','Auto', 'MaxOrder',3, ...
                 'InitialStep','Auto', 'FixedStep','Auto');

EstOpts = nlgreyestOptions('Display','on', ...
                           'SearchMethod','lsqnonlin', ...
                           'EstimateCovariance',false);

Advanced = optimset('lsqnonlin');
% Advanced.Algorithm = 'trust-region-reflective';
Advanced.TolFun = 1e-5;
Advanced.TolX = 1e-6;
Advanced.MaxIter = 100;
EstOpts.GradientOptions.MinDifference = 1e-3;
EstOpts.SearchOptions.Advanced = Advanced;
EstOpts.SearchOptions.FunctionTolerance = Advanced.TolFun;
EstOpts.SearchOptions.StepTolerance = Advanced.TolX;
EstOpts.SearchOptions.MaxIterations = Advanced.MaxIter;
% EstOpts.Regularization.Nominal = 'model';
% EstOpts.Regularization.Lambda = 0.01; % Not directly used by nlgreyest but for reference
% EstOpts.Regularization.R = 1e-3; % For sensitivity analysis post-estimation
% EstOpts.Regularization.R = 1e-3; % For sensitivity analysis post-estimation

EstOpts.GradientOptions.MinDifference = 1e-5;

end
