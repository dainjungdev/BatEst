function [V, S, dSOC_dV, real_dSOC_dV] = load_dSOC_dV(OCP_filename)
% dQ/dV
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

%% 2. Define Voltage function(x:Voltage, y:SOC)
% Plotting
if generate_plot
    figure;
    plot(V, S);
    xlabel('Voltage (V)');
    ylabel('State of charge (SOC)');
    title('State of Charge vs. Voltage');
end

%% 2-1. Downsample voltage
lt = 2000;  % about 2000 data points
ds = max(floor(length(V)/lt),1);  % Example downsampling factor, adjust as needed
downsampledIndices = 1:ds:length(V); % Indices of the downsampled data

% Downsampled data
V = V(downsampledIndices);
S = S(downsampledIndices);

% Make sure that 'V' has unique value
% Aggregate S values for each unique V
[V, ~, ic] = unique(V); % Find unique V values and their indices
S = accumarray(ic, S, [], @mean); % Calculate average S for each unique V

% Define the limits
V_min = 2.5;
V_max = 4.18;

% Find indices where V is within the specified range
valid_indices = V >= V_min & V <= V_max;

% Limit V and S to these indices
V = V(valid_indices);
S = S(valid_indices);

V_S = @(v) ppval(spline(V, S), v);
V = linspace(min(V), max(V), 1000);
S = V_S(V);

% %% Downsample voltage
% lt = 1500;  % about 1500 data points
% ds = max(floor(length(V)/lt),1);  % Example downsampling factor, adjust as needed
% % S = movmean(S, ds);
% downsampledIndices = 1:ds:length(V); % Indices of the downsampled data
% 
% % Downsampled data
% V = V(downsampledIndices);
% S = S(downsampledIndices);

%% 3. Calculate the gradient dSOC_dV
dSOC_dV = gradient(S, V);
real_dSOC_dV = dSOC_dV;

% Plotting
if generate_plot
    figure;
    plot(V, dSOC_dV);
    xlabel('Voltage (V)');
    ylabel('dSOC/dV');
    title('dSOC/dV vs. Voltage');
end

%% 3-1. Gaussian Filter
sigma = 10 + 1; % Standard deviation for Gaussian kernel
dSOC_dV = imgaussfilt(dSOC_dV, sigma);

% Plotting
if generate_plot
    figure;
    TF = islocalmin(dSOC_dV, 'MaxNumExtrema', 4);
    plot(V, dSOC_dV, V(TF), dSOC_dV(TF), 'r*');
    xlabel('Voltage (V)');
    ylabel('dSOC/dV');
    title('dSOC/dV(Gaussian Filtered) vs. Voltage');
end
end