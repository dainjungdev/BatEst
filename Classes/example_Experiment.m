% Create an experiment
exp = Experiment();

% Add instructions including a sine wave
exp = exp.addCurrent(1, 3600);  % Current of 1 A for 3600 seconds
exp = exp.addRest(1800);        % Rest for 1800 seconds
exp = exp.addSineWave(0.5, 1/60, 0, 3600); % Sine wave with amplitude 0.5 A, frequency 1/60 Hz for 3600 seconds
exp = exp.addVoltage(4.1, 1200);% Voltage of 4.1 V for 1200 seconds

% Run the experiment
[out, params] = exp.run();
