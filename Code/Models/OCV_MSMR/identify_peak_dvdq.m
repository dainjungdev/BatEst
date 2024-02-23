function [V, S, dV_dSOC] = identify_peak_dvdq(OCP_filename)
% Identify peaks

reset_path;
close all;
generate_plot = true;

%% 1. Load S & V
% Load the dataset
T = parquetread(OCP_filename{1});
% Ensure that all values are of numeric type double
T = convertvars(T,T.Properties.VariableNames,'double');

% Define starting stoichiometry value
if mean(T.Current_A)<0
    % Discharge
    Start = 0;
else
    % Charge
    Start = 1;
end

% Compute the capacity and stoichiometry/SOC by integration
S = cumtrapz(T.Test_Time_s,T.Current_A);
S = Start+(-1)^Start*S/S(end);
V = T.Voltage_V;

%% 2. Find unique V values to generate plots
% Make sure that 'V' has unique value
% Aggregate S values for each unique V
[V, ~, ic] = unique(V); % Find unique V values and their indices
S = accumarray(ic, S, [], @mean); % Calculate average S for each unique V

%% 3. Filter V and S, to consider voltage within the range
voltage_range = [3.25, 4.1];
indicesInRange = V >= voltage_range(1) & V <= voltage_range(2);

% Apply the filter to get V and dV_dSOC within the desired range
V = V(indicesInRange);
S = S(indicesInRange);

%% 4. Downsample voltage
% Step 1: Downsampling
lt = 1500;  % about 1500 data points
ds = max(floor(length(V)/lt),1);  % Example downsampling factor, adjust as needed
downsampledIndices = 1:ds:length(V); % Indices of the downsampled data

% Downsampled data
V = V(downsampledIndices);
S = S(downsampledIndices);

%% 5. Calculate the gradient dV_dSOC
dV_dSOC = gradient(V, S);

% Now, plotting dV_dSOC vs. Voltage only for the specified range
if generate_plot
    figure;
    plot(V, dV_dSOC);
    xlabel('Voltage (V)');
    ylabel('dV/dSOC');
    xlim(voltage_range);  % Set the x-axis limits to the specified range
    title('dV/dSOC vs. Voltage within [3.25, 4.1]');
end

%% 6. Gaussian Filter
sigma = ds * 2 + 1; % Standard deviation for Gaussian kernel
dV_dSOC = imgaussfilt(dV_dSOC, sigma);

% Plotting
if generate_plot
    figure;
    TF = islocalmax(dV_dSOC, 'MaxNumExtrema', 4);
    plot(V, dV_dSOC, V(TF), dV_dSOC(TF), 'r*');
    xlabel('Voltage (V)');
    ylabel('dSOC/dV');
    xlim(voltage_range);  % Set the x-axis limits to the specified range
    title('dSOC/dV(Gaussian Filtered) vs. Voltage');
end