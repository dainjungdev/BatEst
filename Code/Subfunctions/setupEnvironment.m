function startTime = setupEnvironment()
    close all;
    reset_path;
    startTime = datetime('now');
    fprintf('\nmain_multi started at %s\n', startTime);
end