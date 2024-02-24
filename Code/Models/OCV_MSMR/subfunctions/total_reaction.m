function [x, dx_du] = total_reaction(MSMR_parameters)
num_reactions = size(MSMR_parameters, 1);
T = 25+273.15;

% Initialise arrays for Xj and dXj_dU
xj = cell(num_reactions, 1);
dxj_du = cell(num_reactions, 1);

for i = 1:num_reactions
    % Get parameters for this reaction
    U0 = MSMR_parameters(i, 1);
    Xj_tot = MSMR_parameters(i, 2);
    omega = MSMR_parameters(i, 3);
    
    % Get function handles for Xj and dXj_dU
    [xj{i}, dxj_du{i}] = individual_reactions_function(U0, Xj_tot, omega, T);
end

% Add the function handles
x = @(x) xj{1}(x) + xj{2}(x) + xj{3}(x) + xj{4}(x);
dx_du = @(x) dxj_du{1}(x) + dxj_du{2}(x) + dxj_du{3}(x) + dxj_du{4}(x);
end