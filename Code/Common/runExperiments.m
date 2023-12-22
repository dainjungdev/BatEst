function [out, params] = runExperiment(expObject)
    % Check if the input is an Experiment object
    if ~isa(expObject, 'Experiment')
        error('Input must be an Experiment object.');
    end

    % Extract instructions from the Experiment object
    instructions = expObject.Instructions;

    % Initialize necessary variables
    Dataset = []; % Load or create the dataset as required
    out = [];     % Initialize the output structure
    params = [];  % Load or create initial parameters

    % You may define these or make them inputs to the function
    rep_num = 1:3; % Repetition numbers for the experiment
    outputPath = './BatEst/Data/out/'; % Output path for results

    % Convert instructions to protocol
    if ~isempty(instructions)
        [tt, uu] = convertInstructionsToProtocol(instructions, params);
        params.tt = tt; % Time points
        params.uu = uu; % Input values corresponding to time points
    end

    % Run the main multi-cycle experiment
    [out, params] = main_multi(Dataset, out, params, rep_num, outputPath);
end
