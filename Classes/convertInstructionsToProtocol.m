function [tt, uu] = convertInstructionsToProtocol(instructions, params)
% convert the structured instructions into tt and uu arrays, 
% which represent time points and corresponding input values for the battery experiment. 
% This conversion depends on the nature of each instruction (current, voltage, rest) 
% and the parameters of the battery model.

[mn, hr, Crate, CtoK] = struct2array(params, {'mn','hr','Crate','CtoK'});

% Initialize arrays
tt = [];
uu = [];

% Current time
currentTime = 0;

for i = 1:length(instructions)
    instruction = instructions(i);

    switch instruction.type
        case 'current'
            % Define the time and current amplitude
            newTimes = linspace(currentTime, currentTime + instruction.duration, 100)';
            newCurrents = ones(size(newTimes)) * instruction.amplitude / Crate;

            % Append to arrays
            tt = [tt; newTimes];
            uu = [uu; newCurrents];

        case 'voltage'
            % Define the time and voltage
            newTimes = linspace(currentTime, currentTime + instruction.duration, 100)';
            newVoltages = ones(size(newTimes)) * instruction.voltage / Vmax;

            % Append to arrays
            tt = [tt; newTimes];
            uu = [uu; newVoltages];

        case 'rest'
            % Define the rest period
            newTimes = linspace(currentTime, currentTime + instruction.duration, 100)';
            newCurrents = zeros(size(newTimes)); % Rest implies no current

            % Append to arrays
            tt = [tt; newTimes];
            uu = [uu; newCurrents];

        % Add more cases as needed
    end

    % Update the current time
    currentTime = tt(end);
end

% Assuming current is the only input. If there are more inputs (like temperature),
% additional columns should be added to uu.
end