function params = model_parameters(params,j)
% Set the values of the model parameters and store in the params structure.

% Define model type
Type = 'OCV_MSMR';

% Load cell properties
[hr] = struct2array(params, {'hr'});

% Set model parameters
Qn = 3.3*hr;         % negative electrode capacity (As)
nu = 1;              % negative/positive electrode capacity ratio (non-dim.)
miu = 1;             % cyclable lithium/positive electrode capacity ratio (non-dim.)

% Set constants
alph = 0.5;          % charge transfer coefficients (non-dim.)
Faraday = 96487;     % Faraday's constant (C mol-1)
Rg = 8.314472;       % gas constant (J mol-1 K-1)

% Update capacity
Q = Qn;              % cell capacity (As)


%% Compile all parameters into the params structure
vars = setdiff(who,{'params','vars'});
for i=1:length(vars), params.(vars{i}) = eval(vars{i}); end


end
