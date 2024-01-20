function [tt, uu] = cell_protocol(params)
% A function to manually define or generate a protocol from data. This
% function is called by set_protocol.m. The outputs are the vector of time
% points tt and corresponding input values uu. There may be up to 3 inputs:
% current (A), temperature (deg. C) and measured voltage (V) in that order.

% Load parameters
[mn, hr, Crate, CtoK] = struct2array(params, {'mn','hr','Crate','CtoK'});

% Define the protocol
t_end = 30*mn; % time period (s)
u1 = @(t) 0.5*Crate; % current (A)
u2 = @(t) 25; % temperature (deg. C)
u3 = @(t) 3.5+0.5*t/t_end; % voltage (V)

% Compute the discrete time protocol
nt = 100;
tt(1:nt,1) = linspace(0,t_end,nt)';
uu(1:nt,1) = u1(tt);
uu(1:nt,2) = u2(tt);
uu(1:nt,3) = u3(tt);

% Or, load the protocol from file
%     load('Data/Examples/drive_cycle.parquet','time','current');
%     tt = time;
%     uu = current;


end


% function [tt, uu] = cell_protocol(params)
%     % Load parameters
%     [mn, hr, Crate, CtoK] = struct2array(params, {'mn', 'hr', 'Crate', 'CtoK'});
% 
%     % Define the protocol duration
%     t_end = 30 * mn; % 30 minutes in seconds
% 
%     % Define dynamic charging profile functions
%     frequency = 0.1; % Frequency for sine wave (in Hz)
%     amplitude = 0.5 * Crate; % Amplitude of the sine wave (half of Crate)
%     u1 = @(t) amplitude * sin(2 * pi * frequency * t); % Sine wave for current (A)
%     u2 = @(t) 25; % Constant temperature (deg. C)
%     u3 = @(t) 3.5 + 0.5 * t / t_end; % Linearly increasing voltage (V)
% 
%     % Compute the discrete time protocol
%     nt = 100; % Number of time steps
%     tt = linspace(0, t_end, nt)'; % Time vector
% 
%     % Compute the input values for each time step
%     uu = zeros(nt, 3); % Initialize uu matrix
%     for i = 1:nt
%         uu(i, 1) = u1(tt(i)); % Current
%         uu(i, 2) = u2(tt(i)); % Temperature
%         uu(i, 3) = u3(tt(i)); % Voltage
%     end
% 
%     % Ensure current is non-negative (charging only)
%     uu(:, 1) = max(uu(:, 1), 0);
% end
