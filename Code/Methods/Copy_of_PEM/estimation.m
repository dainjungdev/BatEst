function sol = runGlobalEstimation(Mass, dxdt, yeqn, params)
% Unpack parameters
[X0, tt, U, yn, c0, c_ind, lb, ub, fiX, verbose] = ...
    struct2array(params, {'X0','tt','U','yy','c0','c_ind','lb','ub', ...
                          'fiX','verbose'});

if ~any(yn)
    error('PEM is only able to estimate parameters from data.');
end
if isempty(c_ind)
    error(['At least one parameter must be unknown. Please check the ' ...
           'value of the uncertainties in set_model.m or inform_params.m.']);
end

% Set bounds
lb = lb(c_ind);
ub = ub(c_ind);
lgap = c0(c_ind)-lb;
ugap = ub-c0(c_ind);
flexible_bounds = true;

% Create iddata object to store output, input and sample rate
timestep = tt(2)-tt(1);
data = iddata(yn,[],timestep);
data.Tstart = tt(1);
data.TimeUnit = '';
    
% Define an initial guess for the parameters
for i = 1:length(c_ind)
    parameters(i).Name = ['p' num2str(i)];
    parameters(i).Unit = '';
    parameters(i).Value = c0(c_ind(i));
    parameters(i).Minimum = lb(i);
    parameters(i).Maximum = ub(i);
    parameters(i).Fixed = false;
end

% PEM options
[SimOpts, EstOpts] = DefaultPEMOpts;
if ~verbose
    EstOpts.Display = 'off';
end

% Configure the nonlinear grey-box model
file_name = ['greybox' num2str(length(c_ind))];
order = [size(yn,2) 0 length(X0)];
initial_states = X0;
Ts = 0; % sample time of discrete model or 0 for a continuous time model

% Define the grey-box model
greyModel = idnlgrey(file_name,order,parameters,initial_states,Ts, ...
                    'FileArgument',{dxdt,yeqn},'SimulationOptions',SimOpts);
% greyModel.TimeUnit = '';
% greyModel = setpar(greyModel, 'Minimum', lb, 'Maximum', ub); % Set parameter bounds
% Set whether the initial states are fixed or free estimation parameters
setinit(greyModel,'Fixed',fiX);

% Define optimization options
opt = nlgreyestOptions;
opt.SearchMethod = 'auto'; % Let MATLAB decide the best method
    
% Setup Global Optimization
ms = MultiStart('UseParallel', false, 'StartPointsToRun', 'bounds');
gs = GlobalSearch(ms);
problem = createOptimProblem('fmincon', 'objective', ...
    @(p) computeObjective(p, greyModel, data), 'x0', [parameters.Value], ...
    'lb', lb, 'ub', ub);
    
% Solve the problem
[bestParameters, bestObjective] = run(gs, problem);

% Update model with best parameters found
greyModel = setpar(greyModel, 'Value', num2cell(bestParameters));

% % Pack up solution
% sol.model = finalModel;
% sol.parameters = bestParameters;
% sol.objective = bestObjective;
% 
% % Perform final estimation with nlgreyest using the best parameters as initial guesses
% finalModel = nlgreyest(data, greyModel, opt);

% Estimate the model parameters and initial states
runtic = tic;
sys = pem(data,greyModel,EstOpts);
X0 = sys.Report.Parameters.X0; %findstates(sys,data);
toc(runtic);

% Extract the parameter estimates
ce = sys.Report.Parameters.ParVector;

% Update the parameters
c0(c_ind) = ce;

% Pack up solution
sol.tsol = tt;
sol.xsol = X0';
sol.usol = U(tt);
sol.psol = ones(length(tt),1)*c0';
end

function obj = computeObjective(parameters, model, data)
    % Update model with new parameters
    updatedModel = setpar(model, 'Value', num2cell(parameters));
    sys = pem(data, updatedModel);
    % Compute objective (e.g., sum of squared errors)
    obj = sys.Report.Fit.LossFcn;
end
