function finalizeComputation(startTime, outputPath, ModelName, out)
    endTime = datetime('now');
    duration = endTime - startTime;
    fprintf('\nComputation ended at %s\n', endTime);
    fprintf('Total duration: %s\n', duration);

    save_output(out, [outputPath '/' ModelName], true);
    reset_path;
end