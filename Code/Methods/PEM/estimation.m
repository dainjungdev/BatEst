function sol = estimation(Mass,dxdt,yeqn,params)
% Perform state/parameter estimation using PEM.

% Unpack parameters
[X0, tt, U, yn, c0, c_ind, lb, ub, fiX, verbose, DataType] = ...
    struct2array(params, {'X0','tt','U','yy','c0','c_ind','lb','ub', ...
                          'fiX','verbose','DataType'});

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

% Options for main_multi
if strcmp(DataType,'Pseudo-OCV charge')
    EstOpts.SearchOptions.StepTolerance = 1e-6;
    EstOpts.SearchOptions.FunctionTolerance = 1e-12;
elseif strcmp(DataType,'Relaxation')
elseif strcmp(DataType,'CCCV charge')
    EstOpts.SearchOptions.StepTolerance = 1e-4;
    EstOpts.SearchOptions.FunctionTolerance = 1e-14;
end

if ~verbose
    EstOpts.Display = 'off';
end

% Configure the nonlinear grey-box model
file_name = ['greybox' num2str(length(c_ind))];
order = [size(yn,2) 0 length(X0)]; % number of outputs, inputs and states
initial_states = X0; % warm start
Ts = 0; % sample time of discrete model or 0 for a continuous time model
init_sys = idnlgrey(file_name,order,parameters,initial_states,Ts, ...
                    'FileArgument',{dxdt,yeqn},'SimulationOptions',SimOpts);
init_sys.TimeUnit = '';

% Set whether the initial states are fixed or free estimation parameters
setinit(init_sys,'Fixed',fiX);

for step = 1
    % Estimate the model parameters and initial states
    runtic = tic;
    sys = pem(data,init_sys,EstOpts);

    % Display key optimization results
    fprintf('\n***** PEM Optimization Results *****\n');
    fprintf('* Loss Function Value:      %.6e\n', sys.Report.Fit.LossFcn);
    fprintf('* Termination Reason:       %s\n', sys.Report.Termination.WhyStop);
    fprintf('* First Order Optimality:   %.6e\n', sys.Report.Termination.FirstOrderOptimality);
    fprintf('* Function Evaluation Count: %d\n', sys.Report.Termination.FcnCount);
    fprintf('************************************\n\n');

    X0 = sys.Report.Parameters.X0; %findstates(sys,data);
    toc(runtic);
    
    % Extract the parameter estimates
    ce = sys.Report.Parameters.ParVector;
    % 
    % lb = sqrt(lb .* ce)
    % ub = sqrt(ce .* ub)
    % % 0.05 0.5 5.0 1
    % % 0.1581 0.5 1.581 2
    % % 0.2811 0.5 0.8891 3
    % % 0.4714 0.5 0.6667 4
    % % 0.4854 0.5 0.5773 5
    % % 0.4926 0.5 0.5372 6
    % % 0.4962 0.5 0.5182 7
    % % 0.4980 0.5 0.5090 8
    % % 0.4989 0.5 0.5044 9
    % % 0.4994 0.5 0.5021 10
    % init_sys = setpar(init_sys,'Value',  num2cell(ce));
    % init_sys = setpar(init_sys,'Minimum',num2cell(lb));
    % init_sys = setpar(init_sys,'Maximum',num2cell(ub));
    % 
    % % Update the initial states
    % init_sys = setinit(init_sys,'Value',num2cell(X0));
    
end

% Update the parameters
c0(c_ind) = ce;

% Pack up solution
sol.tsol = tt;
sol.xsol = X0';
sol.usol = U(tt);
sol.psol = ones(length(tt),1)*c0';

end