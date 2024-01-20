function startTime = setupEnvironment()
    close all;
    reset_path;
    startTime = datetime('now');
    fprintf('\nCode started at %s\n', startTime);
end