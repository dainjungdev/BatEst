function [true_sol, params] = step11(Target, Model, params, j, Dataset, experiment)
% This step loads the data or simulates the model to generate data.

% Initialise the solution structure
true_sol = [];

% Unpack model properties
[ModelName, Noise] = struct2array(Model, {'Name', 'Noise'});

%% Check if Dataset is provided and Target is one of the specified types
if istable(Dataset) && any(strcmp(Target,{'Simulate','Plot','Compare','Parameter'}))
    % Unpack the data from the dataset
    true_sol = unpack_data(Dataset, params, j);
    true_sol.Type = 'True';

    % Update the parameters accordingly
    params = inform_params(params, ModelName, true_sol);
    
    % Set the protocol based on the data
    params = set_protocol(params, true_sol);

% Set the protocol based on script if no dataset is provided
elseif any(strcmp(Target,{'Simulate', 'Parameter'}))

    % Handle the case where an Experiment object is provided
    if nargin == 6 && isa(experiment, 'Experiment')
        % Set the protocol based on the instructions from the Experiment object
        params = set_protocol(params, experiment);

    else
        params = set_protocol(params);
    
    end
    
% If Target is not recognised, throw an error
elseif ~strcmp(Target,'Control')
    error(['Check that Target is set to one of the available options ' ...
           'and/ or import a Dataset and pass it to the function by ' ...
           'entering main(Dataset); in the command window.']);
end


%% If simulation is required or no dataset is provided for parameter setting
if strcmp(Target,'Simulate') || (~istable(Dataset) && strcmp(Target,'Parameter'))
    % Run the model simulation
    disp('hi');
    true_sol = run_simulation(Model, params);
    true_sol.Type = 'True';  % Mark the solution as true data
end


%% If a solution structure is available, plot and save the results
if isstruct(true_sol)
    % Plot any simulation
    params = plot_sol(true_sol, params);
    % Option to add measurement noise with a given standard deviation
    if Noise
        std = 0.01;
        true_sol.ysol = true_sol.ysol+std*randn(size(true_sol.ysol));
        params = plot_sol(true_sol,params);
    end

    % Save the solution data
    params.yy = true_sol.ysol;
    if params.verbose
        disp('Forward pass complete.');
    end
end

end
