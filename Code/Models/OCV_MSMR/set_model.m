function [Model, params] = set_model(ModelName,params,j)
% This function defines an equivalent circuit model with an OCV-MSMR source.
% A positive current (I>0) corresponds to charging.
% Need

% The model is given by:
% dSOC/dt = I/Q,
% V = OCV(SOC).

% The model is scaled and written below in terms of:
% dimensionless time t/Tm,
% the input  u = [I/Um],
% the states x = [SOC],
% the output y = [(V-Vcut)/Vrng].

% Unpack parameters
[Q, nu, miu, Um, Vrng, Vcut, OCV, Tm, S0] = ...
    struct2array(params, {'Q','nu','miu',...
                          'Um','Vrng','Vcut','OCV','Tm','S0'});

% Define an initial guess and uncertainty for each unknown parameter
guess = [1/Q; nu; miu];
uncert = [0.1; 0.1; 0.1];

% Set the rescaling factor and scale the initial guesses
fac = 2*guess;
c0 = guess./fac;

% write OCV as a function of c
% number of equations: number of peaks appearing in the derivatives



% Compile parameters into vector
c = {c0(1); c0(2); c0(3); ... [1-3] length of guess
     Um; Vcut; Vrng; ... scaling [4-6]
     OCV; ... subfunctions [7]
     Tm}; % keep the timescale Tm as the last entry [8]

% Define the number of parameters (not including subfunctions or Tm)
params.nop = 6;

% Define helper function
f = @(c,i,t) feval(c{i},t)*fac(i);

% Define the state derivatives
dxdt = @(t,x,y,u,c) [c{1}*c{4}*u(1,:)]*c{8};

% Define the output equation
yeqn = @(t,x,u,c) (c{7}(x(1,:),c{2},c{3})-c{5})/c{6};

% Define the mass matrix
Mass = diag([1; 0]);

% Set the initial states
params.X0 = [S0];

% Pack up the model
Model = struct('Name', ModelName, 'Mass',Mass, 'dxdt',dxdt, 'yeqn', yeqn);
params.uncert = uncert; params.fac = fac; params.c0 = c0; params.c = c;

end
