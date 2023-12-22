classdef Experiment
    properties
        Instructions % Array to store the instructions
    end

    methods
        function obj = Experiment()
            % Constructor for the Experiment class
            obj.Instructions = {};
        end

        function obj = addCurrent(obj, amplitude, duration)
            % Add a current instruction
            obj.Instructions{end+1} = current(amplitude, duration);
        end

        function obj = addVoltage(obj, volt, duration)
            % Add a voltage instruction
            obj.Instructions{end+1} = voltage(volt, duration);
        end

        function obj = addRest(obj, duration)
            % Add a rest instruction
            obj.Instructions{end+1} = rest(duration);
        end
        
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



        end
    end
end
