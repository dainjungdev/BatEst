function params = model_parameters(params,j)
% Set the values of the model parameters and store in the params structure.

% Define model type
Type = 'EHM';

% Load cell properties
[hr, Tamb, Tref] = struct2array(params, {'hr','Tamb','Tref'});

% Set model parameters
Qn = 3.3*hr;         % negative electrode capacity (As)
nu = 1;              % negative/positive electrode capacity ratio (non-dim.)
miu = 1;             % cyclable lithium/positive electrode capacity ratio (non-dim.)
tau_ref = 200;       % diffusion time constant (s) 200
b = 0.3;             % negative electrode surface/particle volume ratio (non-dim.)
Ip_ref = 20;         % reference exchange current in the positive electrode (A)
In_ref = 2;          % reference exchange current in the negative electrode (A) 2
Rf = 0.01;           % film resistance (Ohm) 0.01

% Update capacity
CE = 1;              % coulombic efficiency (non-dim.)
Q = Qn*CE;           % effective negative electrode capacity (As)

% Set constants
alph = 0.5;          % charge transfer coefficients (non-dim.)
Faraday = 96487;     % Faraday's constant (C mol-1)
Rg = 8.314472;       % gas constant (J mol-1 K-1)


%% Temperature dependence
% Activation energies (J mol-1)
E_Dsn = 42770; % for solid-state diffusion in the negative electrode
E_kn  = 37480; % for the reaction in the negative electrode
E_kp  = 39570; % for the reaction in the positive electrode

% Activated parameter values
tau = tau_ref/exp(E_Dsn/Rg*(1/Tref-1/Tamb)); % diffusion time constant (s)
Ip = Ip_ref*exp(E_kp/Rg*(1/Tref-1/Tamb)); % maximum exchange current (A)
In = In_ref*exp(E_kn/Rg*(1/Tref-1/Tamb)); % maximum exchange current (A)


%% Compile all parameters into the params structure
vars = setdiff(who,{'params','vars'});
for i=1:length(vars), params.(vars{i}) = eval(vars{i}); end


end

