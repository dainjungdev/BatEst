function [V, S, dSOC_dV, real_dSOC_dV, V_ext, dSOC_dV_ext] = load_dSOC_dV_edit(OCP_filename)
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

%% 2. Define Voltage function
% Plotting (x:Voltage, y:SOC)
if generate_plot
    figure;
    plot(V, S);
    xlabel('Voltage (V)');
    ylabel('State of charge (SOC)');
    title('State of Charge vs. Voltage');
end

% Plotting (x:SOC, y:Voltage)
if generate_plot
    figure;
    plot(S, V);
    xlabel('State of charge (SOC)');
    ylabel('Voltage (V)');
    title('Voltage vs. State of Charge');
end

%% 3. Pre-processing
%% Filter non-unique indices
% Find unique values in V and their indices, and count of each V
[S_unique, ~, ic] = unique(S);
counts = accumarray(ic, 1);

% Find V values that occur exactly once
S_single_occurrence = S_unique(counts == 1);

% Find indices of V that correspond to V values occurring exactly once
idx_single_occurrence = ismember(S, S_single_occurrence);

% Filter V and S to keep only those with V values occurring exactly once
V = V(idx_single_occurrence);
S = S(idx_single_occurrence);

%% Downsampling
spl_s_v = pchip(S, V);
s_v = @(S) ppval(spl_s_v, S);
% S = linspace(0, 1, 5000);
% V = s_v(S);

ds = ceil(length(S)/5000);

S = [S(1:ds:end-1); S(end)];
V = [V(1:ds:end-1); V(end)];

%% Filter non-monotonic indices
for j = 1:5
% Initialize an empty array to store indices
nonMonotonicIndices = [];

% Loop through the vector starting from the second element
for i = 2:length(V)-1  
    if (V(i-1)<V(i) && V(i+1)<V(i) && V(i-1)>V(i+1)) || (V(i-1)>V(i) && V(i+1)>V(i) && V(i-1)>V(i+1)) 
        nonMonotonicIndices = [nonMonotonicIndices, i]; % Store the index
    end
end

if isempty(nonMonotonicIndices) break; end

% Create a logical index array with all elements set to true initially
logicalIndex = true(1, length(V));

% Set elements at nonMonotonicIndices to false
logicalIndex(nonMonotonicIndices) = false;

% Use the logical index to retain only the elements that are true (not in nonMonotonicIndices)
V = V(logicalIndex);
S = S(logicalIndex);
end

%% Filter non-unique indices
% Find unique values in V and their indices, and count of each V
[V_unique, ~, ic] = unique(V);
counts = accumarray(ic, 1);

% Find V values that occur exactly once
V_single_occurrence = V_unique(counts == 1);

% Find indices of V that correspond to V values occurring exactly once
idx_single_occurrence = ismember(V, V_single_occurrence);

% Filter V and S to keep only those with V values occurring exactly once
V = V(idx_single_occurrence);
S = S(idx_single_occurrence);

%% Filter non-monotonic indices
for j = 1:5
% Initialize an empty array to store indices
nonMonotonicIndices = [];

% Loop through the vector starting from the second element
for i = 2:length(V)-1  
    if (V(i-1)<V(i) && V(i+1)<V(i) && V(i-1)>V(i+1)) || (V(i-1)>V(i) && V(i+1)>V(i) && V(i-1)>V(i+1)) 
        nonMonotonicIndices = [nonMonotonicIndices, i]; % Store the index
    end
end

if isempty(nonMonotonicIndices) break; end

% Create a logical index array with all elements set to true initially
logicalIndex = true(1, length(V));

% Set elements at nonMonotonicIndices to false
logicalIndex(nonMonotonicIndices) = false;

% Use the logical index to retain only the elements that are true (not in nonMonotonicIndices)
V = V(logicalIndex);
S = S(logicalIndex);
end

%% Check limits
% Define the limits
V_min = 2.5;
V_max = 4.18;

% Find indices where V is within the specified range
valid_indices = V >= V_min & V <= V_max;

% Limit V and S to these indices
V = V(valid_indices);
S = S(valid_indices);

V_S = @(v) ppval(pchip(V, S), v);
V = linspace(min(V), max(V), 20001);
S = V_S(V);

% if generate_plot
%     figure;
%     plot(V, S, 'b'); hold on;
%     plot(V, V_S(V),'r'); hold off;
%     title('After downsampling');
% end

%% 3. Calculate Gradient
% Preallocate the gradient array
dSOC_dV_arr = zeros(200, size(S,2));

for j = 1:200
    % Interval size
    interval = j; % Adjust based on your data and desired smoothing
    
    dSOC_dV_arr(j,1) = (S(2) - S(1)) / (V(2) - V(1));
    for i = 2:interval
        dSOC_dV_arr(j,i) = (S(i*2-1) - S(1)) / (V(i*2-1) - V(1));
    end
    
    % Calculate the gradient over the interval
    for i = 1 + interval:length(S) - interval
        dSOC_dV_arr(j,i) = (S(i+interval) - S(i-interval)) / (V(i+interval) - V(i-interval));
    end
    
    for i = length(S)-interval+1:length(S)
        dSOC_dV_arr(j,i) = (S(end) - S(i*2-end)) / (V(end) - V(i*2-end));
    end
    dSOC_dV_arr(j,end) = (S(end) - S(end-1)) / (V(end) - V(end-1));
end

% for i = 1:length(S)
%     dSOC_dV_arr(:, i) = imgaussfilt(dSOC_dV_arr(:, i), 3);
% end

real_dSOC_dV_weights = linspace(1, 0, length(20:200)); % Linear weights from 1 to 0
real_dSOC_dV_weights = real_dSOC_dV_weights / sum(real_dSOC_dV_weights); % Normalize weights
real_dSOC_dV = sum(dSOC_dV_arr(20:200,:).* real_dSOC_dV_weights');

% dSOC_dV_weights = linspace(1, 1.5, length(10:300)); % Linear weights from 1 to 1.5
% dSOC_dV_weights = dSOC_dV_weights / sum(dSOC_dV_weights); % Normalize weights
% dSOC_dV = sum(dSOC_dV_arr(10:300,:).* dSOC_dV_weights');

real_dSOC_dV = mean(dSOC_dV_arr(20:200,:));
dSOC_dV = mean(dSOC_dV_arr(20:200,:));

% dSOC_dV_weights = linspace(2, 1, length(20:200)); % Linear weights from 1 to 1.5
% dSOC_dV_weights = dSOC_dV_weights / sum(dSOC_dV_weights); % Normalize weights
% dSOC_dV = sum(dSOC_dV_arr(20:200,:).* dSOC_dV_weights');


V = V(1:10:end);
S = S(1:10:end);

% 3-1. Gaussian Filter
sigma = 10; % Standard deviation for Gaussian kernel
dSOC_dV = imgaussfilt(dSOC_dV, sigma, 'FilterSize', 21);

dSOC_dV = dSOC_dV(1:10:end);
real_dSOC_dV = real_dSOC_dV(1:10:end);

% Sgolay filter
polyOrder = 5; windowSize = 101;
dSOC_dV = sgolayfilt(dSOC_dV, polyOrder, windowSize);

% Plot the gradient
if generate_plot
    figure;
    TF = islocalmin(dSOC_dV, 'MaxNumExtrema', 4);
    plot(V, real_dSOC_dV, 'r'); hold on;
    plot(V, dSOC_dV, 'g');
    plot(V(TF), dSOC_dV(TF), 'r*');
    xlabel('Voltage (V)');
    ylabel('dSOC/dV');
    title('Calculated Gradient with Increased Interval');
end

% %% 3-1. Gaussian Filter
% sigma = 1; % Standard deviation for Gaussian kernel
% dSOC_dV = imgaussfilt(dSOC_dV, sigma);

% polyOrder = 5; windowSize = 51;
% dSOC_dV = sgolayfilt(dSOC_dV, polyOrder, windowSize);

% % Plotting
% if generate_plot
%     figure;
%     TF = islocalmin(dSOC_dV, 'MaxNumExtrema', 4);
%     plot(V, dSOC_dV, V(TF), dSOC_dV(TF), 'r*'); hold on;
%     plot(V, real_dSOC_dV);
%     xlabel('Voltage (V)');
%     ylabel('dSOC/dV');
%     title('dSOC/dV(Gaussian Filtered) vs. Voltage');
% end

end