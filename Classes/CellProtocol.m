classdef CellProtocol
    properties
        params % Parameters of the cell
    end

    methods
        function obj = CellProtocol(params)
            obj.params = params;
        end

        function [tt, uu] = defaultProtocol(obj)
            % Default protocol generation
            % This method generates a default protocol. This could be a standard
            % charging and discharging cycle or any other protocol you commonly 
            % use.

            t_end = 30 * obj.params.mn; % time period in seconds
            nt = 100; % Number of time points
        
            % Time vector
            tt = linspace(0, t_end, nt)';
        
            % Current (A), Temperature (deg. C), Voltage (V)
            uu = zeros(nt, 3); 
            uu(:, 1) = 0.5 * obj.params.Crate; % current
            uu(:, 2) = 25; % temperature
            uu(:, 3) = 3.5 + 0.5 * (tt / t_end); % voltage
        end

        function [tt, uu] = fromInstructions(obj, instructions)
            % Define protocol from a set of instructions
            % This method takes a set of instructions (similar to those in your Experiment
            % class) and converts them into a protocol with time points tt and corresponding
            % input values uu.
            
            disp('cellprotocoltest')
            disp('hi')
            % Convert structured instructions into tt and uu arrays
            Crate = obj.params.Crate;
            Vmax = obj.params.Vmax; % Assuming these exist in your params structure
        
            % Initialize arrays
            tt = [];
            uu = [];
        
            % Current time
            currentTime = 0;
        
            for i = 1:length(instructions)
                instruction = instructions(i);
                instruction = instruction{1,1};
        
                switch instruction.type
                    case 'current'
                        % Time and current amplitude
                        newTimes = linspace(currentTime, currentTime + instruction.duration, 100)';
                        newCurrents = ones(size(newTimes)) * instruction.amplitude / Crate;
        
                        % Append to arrays
                        tt = [tt; newTimes];
                        uu = [uu; newCurrents, zeros(size(newTimes)), zeros(size(newTimes))];
        
                    case 'voltage'
                        % Time and voltage
                        newTimes = linspace(currentTime, currentTime + instruction.duration, 100)';
                        newVoltages = ones(size(newTimes)) * instruction.voltage / Vmax;
        
                        % Append to arrays
                        tt = [tt; newTimes];
                        uu = [uu; zeros(size(newTimes)), zeros(size(newTimes)), newVoltages];
        
                    case 'rest'
                        % Rest period
                        newTimes = linspace(currentTime, currentTime + instruction.duration, 100)';
                        newCurrents = zeros(size(newTimes)); % Rest implies no current
        
                        % Append to arrays
                        tt = [tt; newTimes];
                        uu = [uu; newCurrents, zeros(size(newTimes)), zeros(size(newTimes))];

                    case 'sineWave'
                        % Define the time and current for sine wave
                        newTimes = linspace(currentTime, currentTime + instruction.duration, 100)';
                        newCurrents = instruction.amplitude * sin(2 * pi * instruction.frequency * newTimes + instruction.phase);
                
                        % Append to arrays
                        tt = [tt; newTimes];
                        uu = [uu; newCurrents];
                
                        % Update the current time
                        currentTime = tt(end);

                    % Add more cases as needed
                end
        
                % Update the current time
                currentTime = tt(end);
            end
        
            % Assuming current is the only input. If there are more inputs (like temperature),
            % additional columns should be added to uu.
        end

        % Additional methods for other specific protocols can be added here
    end
end
