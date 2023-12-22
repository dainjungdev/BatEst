function instruction = current(amplitude, duration)
    % Create a current instruction
    instruction = struct('type', 'current', 'amplitude', amplitude, 'duration', duration);
end