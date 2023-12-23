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
            obj.Instructions{end+1} = struct('type', 'current', 'amplitude', amplitude, 'duration', duration);
        end

        function obj = addVoltage(obj, volt, duration)
            % Add a voltage instruction
            obj.Instructions{end+1} = struct('type', 'voltage', 'voltage', volt, 'duration', duration);
        end

        function obj = addRest(obj, duration)
            % Add a rest instruction
            obj.Instructions{end+1} = struct('type', 'rest', 'duration', duration);
        end

        function obj = addSineWave(obj, amplitude, frequency, phase, duration)
            % Add a sine wave instruction
            sineWave = struct('type', 'sineWave', 'amplitude', amplitude, 'frequency', frequency, 'phase', phase, 'duration', duration);
            obj.Instructions{end+1} = sineWave;
        end

        function [out, params] = run(obj, params, outputPath)
            % Run the experiment
            if isempty(obj.Instructions)
                error('No instructions have been added to the experiment.');
            end

            % Set default values for optional parameters
            if nargin < 2, params = []; end
            if nargin < 3, outputPath = './BatEst/Data/out/'; end

            % Define repetition numbers and dataset
            rep_num = 1:3;
            Dataset = []; % Or load/create dataset as required

            % Run the main multi-cycle experiment
            [out, params] = main_multi_experiment(Dataset, [], params, rep_num, outputPath, obj);
        end
    end
end
