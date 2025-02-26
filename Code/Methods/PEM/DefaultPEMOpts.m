function [SimOpts, EstOpts] = DefaultPEMOpts
% Set the default simulation options for PEM.

SimOpts = struct('Solver','Auto', 'RelTol',1e-3, 'AbsTol',1e-5, ...
                 'MinStep','Auto', 'MaxStep','Auto', 'MaxOrder',5, ...
                 'InitialStep','Auto', 'FixedStep','Auto');

EstOpts = nlgreyestOptions('Display','on', ...
                           'SearchMethod','lsqnonlin', ...
                           'EstimateCovariance',true);

Advanced = optimset('lsqnonlin');
% Advanced.Algorithm = 'trust-region-reflective';

Advanced.TolFun = 1e-12;
Advanced.TolX = 1e-4;
Advanced.MaxIter = 100;

EstOpts.SearchOptions.Advanced = Advanced;
EstOpts.EstimateCovariance = true;
EstOpts.SearchOptions.FunctionTolerance = Advanced.TolFun;
EstOpts.SearchOptions.StepTolerance = Advanced.TolX;
EstOpts.SearchOptions.MaxIterations = Advanced.MaxIter;

EstOpts.GradientOptions.MinDifference = 1e-4;

end