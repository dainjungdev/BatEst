% A script to test runExperiment
% Create an experiment with a mix of current, voltage, and rest instructions
instructions = [
    current(1, 3600),  % Discharge at 1 A for 1 hour
    rest(1800),        % Rest for 30 minutes
    voltage(4.1, 1800) % Hold at 4.1 V for 30 minutes
];

% Run the experiment
[out, params] = runExperiments(instructions);

% Inspect the output and params for results
disp(out);
disp(params);
