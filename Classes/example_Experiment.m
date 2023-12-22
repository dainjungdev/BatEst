% Create an experiment
exp = Experiment();

% Add instructions
exp = exp.addCurrent(1, 3600);  % Current of 1 A for 3600 seconds
exp = exp.addRest(1800);        % Rest for 1800 seconds
exp = exp.addVoltage(4.1, 1200);% Voltage of 4.1 V for 1200 seconds

% Run the experiment
exp.run();