function [est_sol,params] = step2(Target, Model, params, j)
% This step runs the estimation to obtain control or parameter estimates.

% Initialise the estimated solution structure
est_sol = [];

%% Handle different Target options
    switch Target
        case {'Parameter', 'Control'}
            % If Target is 'Parameter', define bounds for parameter estimation
            if strcmp(Target, 'Parameter')
                params = add_bounds(params);
            end

            % Set up and run the parameter estimation
            [Mass, est_dxdt, est_yeqn, params] = set_unknown(Target, Model, params);
            est_sol = estimation(Mass, est_dxdt, est_yeqn, params);
            est_sol.Type = Target;  % Mark the type of estimation

            % If Target is 'Parameter', update parameters with estimated values
            if strcmp(Target, 'Parameter')
                params = retrieve_params(params, est_sol);
            end

            % Plot any estimated solution
            if isfield(est_sol,'ysol')
                params = plot_sol(est_sol, params);
            end
            if params.verbose
                disp('Estimation complete.');
            end
            
        case {'Simulate', 'Plot', 'Compare'}
            % Nothing to estimate
            
        otherwise
            error('Invalid Target specified for step2 function.');
    end
end
