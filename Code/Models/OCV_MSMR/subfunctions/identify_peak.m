function [V, dSOC_dV, dSOC_dV_V] = identify_peak(OCP_filename)
reset_path;
close all;
generate_plot = false;

%% 1. Define OCV function(x:SOC, y:Voltage)
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

[S, ~, ic] = unique(S); % Find unique S values and their indices
V = accumarray(ic, V, [], @mean); % Calculate average V for each unique S

% Define the OCV function
spl = spline(S, V); % produces a piecewise polynomial for use by PPVAL - can use spline, pchip or makima
OCV = @(S) ppval(spl,S);

% Plot
if generate_plot
    figure;
    subplot(1, 2, 1);
    plot(S, OCV(S));
    xlabel('State of charge (SOC)');
    ylabel('Voltage (V)');
    title('Voltage vs. State of Charge');
end

%% 2. Define Voltage function(x:Voltage, y:SOC)
% Make sure that 'V' has unique value
% 1) Aggregate S values for each unique V
[V, ~, ic] = unique(V); % Find unique V values and their indices
S = accumarray(ic, S, [], @mean); % Calculate average S for each unique V

% 2) Create the spline with processed data
spl = spline(V, S);

% 3) Generate the inverse function
SOC_V = @(V) ppval(spl, V);

% Plotting
if generate_plot
    subplot(1, 2, 2);
    SOC = SOC_V(V);
    plot(V, SOC);
    xlabel('Voltage (V)');
    ylabel('State of charge (SOC)');
    title('State of Charge vs. Voltage');
end

%% 2-1. Downsampling SOC_V
% Step 1: Downsampling
lt = 1500;  % about 1500 data points
ds = max(floor(length(V)/lt),1);  % Example downsampling factor, adjust as needed
S = movmean(S, ds);
downsampledIndices = 1:ds:length(V); % Indices of the downsampled data

% Downsampled data
V = V(downsampledIndices);
S = S(downsampledIndices);

% Step 2: Create the Spline with Downsampled Data
spl = spline(V, S);

% Step 3: Generate the Inverse Function with Downsampled Data
SOC_V = @(V) ppval(spl, V);

% Plotting
figure;
SOC = SOC_V(V);
plot(V, SOC);
xlabel('Voltage (V)');
ylabel('State of charge (SOC)');
title('Downsampled Average State of Charge vs. Voltage');

%% 3. Calculate the gradient dSOC_dV
dSOC_dV = gradient(S, V);

% Plotting
figure;
plot(V, dSOC_dV);
xlabel('Voltage (V)');
ylabel('dSOC/dV');
title('dSOC/dV vs. Voltage');

%% 3-1. Gaussian Filter
sigma = ds * 2 + 1; % Standard deviation for Gaussian kernel
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

%% 3-2. Generate dSOC_dV function
spl = spline(V, dSOC_dV);
dSOC_dV_V = @(V) ppval(spl, V);

% % Save plot
% outputPath = 'Project/07-MSMR/pseudoOCV_analysis/out';
% % mkdir(outputPath);
% fileName = sprintf('Cell%d_OCV%d', cell_num, ocv_num);
% save_plot(gcf, [outputPath '/' fileName], true);

end