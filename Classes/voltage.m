function instruction = voltage(volt, duration)
    % Create a voltage instruction
    instruction = struct('type', 'voltage', 'voltage', volt, 'duration', duration);
end
