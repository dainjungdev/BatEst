function SOC = initial_SOC(params,V,SOCest)
% A function to estimate the steady-state SOC from a known voltage. The
% inputs are the parameters structure params, the voltage V, and the SOC
% estimate SOCest. UpFun and UnFun are the positive and negative half-cell
% OCPs as functions of the SOC (based on the negative electrode).

% Unpack parameters
[nu, miu, UpFun, UnFun, OCV] = ...
    struct2array(params, {'nu','miu','UpFun','UnFun','OCV'});

% Define the function to solve
if isa(OCV,'function_handle')
    eqn = @(x) OCV(x,nu,miu) - V;
else
    eqn = @(x) UpFun(x,nu,miu) - UnFun(x) - V;
end

% Estimate the equilibrium SOC corresponding to the given voltage
opts = optimoptions('fsolve', 'Algorithm', 'trust-region-dogleg', 'Display', 'off');
[SOC, ~, exitflag] = fsolve(eqn, SOCest, opts);

% Retry with 'levenberg-marquardt' if needed
if exitflag <= 0
    opts.Algorithm = 'levenberg-marquardt';
    SOC = fsolve(eqn, SOCest, opts);
end
end