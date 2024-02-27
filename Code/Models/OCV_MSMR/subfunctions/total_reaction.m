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

% Initialize the function handles to zero functions
x = @(x) 0;
dx_du = @(x) 0;

% Loop through each reaction and add the function handles
for i = 1:num_reactions
    % Capture the current state of x and dx_du in the loop to avoid referencing issues
    current_x = x;
    current_dx_du = dx_du;
    x = @(x) current_x(x) + xj{i}(x);
    dx_du = @(x) current_dx_du(x) + dxj_du{i}(x);
end

end