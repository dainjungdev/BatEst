function [MSMR_parameters, RMSE] = get_MSMR_parameters(OCP_filename)
generate_plot = false;

%% 0. Settings
% Define essential constants
T = 298.15; % temp
Rg = 8.314472;  % gas constant (J mol-1 K-1)
F = 96487;     % Faraday's constant (C mol-1)
f = F / (Rg * T);  % F/RT

% Get V, dSOC_dV data
[V, S, dSOC_dV, real_dSOC_dV] = identify_peak(OCP_filename);


%% 1. Get initial estimates
% Initialise the parameter vectors
Q_initial = zeros(1, 4);
U0_initial = zeros(1, 4);
w_initial = zeros(1, 4);

dSOC_dV_to_fit = dSOC_dV;

% Loop over the number of sigmoids to fit
for j = 1:4
    %% 1. Get the first peak
    % Add: Maybe later we can define limits for first peak
    % Invert the signal to find troughs as peaks
    [pks, locs, w, p] = findpeaks(-dSOC_dV_to_fit, 'SortStr', 'descend', 'NPeaks', 5);

    peakLocs = [locs', p'.*w'.*pks', p', w', pks'];
    peakLocs = sortrows(peakLocs, 2, 'descend');
    peakLoc = peakLocs(1);

    V_peak = V(peakLoc);
    dSOC_dV_peak = dSOC_dV_to_fit(peakLoc);

    % Use the location of the first peak as U0
    U0_initial(j) = V_peak;
    

    %% 3. Compute w_initial
    p = linspace(0.95, 0.05, 200); range = norminv(p, 0.5, 0.3);
    p = linspace(0.8, 0.2, 200); range = norminv(p, 0.5, 0.5);

    U0_adjustment = U0_initial(j);
    counter = 0;

    for i = 1:length(range)
        k = range(i);
        ref_dSOC_dV = dSOC_dV_peak * k;
        V_ref_left = findCrossing(V, dSOC_dV_to_fit, peakLoc, ref_dSOC_dV, 'left');
        V_ref_right = findCrossing(V, dSOC_dV_to_fit, peakLoc, ref_dSOC_dV, 'right');

        if ~isnan(V_ref_left) && ~isnan(V_ref_right)
            ref_distance = (V_ref_right - V_ref_left) * 1/2;
            w = ref_distance * f / log(2/k - 1 + 2*sqrt(1/(k*k) - 1/k));
            w_initial(j) = (w_initial(j) * (i - 1) + w) / i;
    
            ref_halfpoint = (V_ref_right + V_ref_left) * 1/2;
            U0_adjustment = (U0_adjustment * (i - 1) + ref_halfpoint) / i;

        elseif ~isnan(V_ref_left)
            ref_distance = V_peak - V_ref_left;
            w = ref_distance * f / log(2/k - 1 + 2*sqrt(1/(k*k) - 1/k));
            counter = counter + 1;
            w_initial(j) = (w_initial(j) * (i - 1/2 - 1/2 * counter) + w * 1/2) / (i - 1/2 * counter);
            if counter >= min(length(range) - i, length(range) / 4) break; end

        elseif ~isnan(V_ref_right)
            ref_distance = V_ref_right - V_peak;
            w = ref_distance * f / log(2/k - 1 + 2*sqrt(1/(k*k) - 1/k));
            counter = counter + 1;
            w_initial(j) = (w_initial(j) * (i - 1/2 - 1/2 * counter) + w * 1/2) / (i - 1/2 * counter);
            if counter >= min(length(range) - i, length(range) / 4) break; end

        else
            break;
        end
    end

    U0_initial(j) = (U0_adjustment + U0_initial(j)) * 1/2;

    % 5. compute Q_initial
    Q_initial(j) = dSOC_dV_peak * (-4) * w_initial(j) / f;
 
    % 6. Get estimated sigmoid
    [Xj, dXj_dU] = individual_reactions(V, U0_initial(j), Q_initial(j), w_initial(j), T); 

    if generate_plot
    figure;
    plot(V, dXj_dU, '--');
    hold on;
    plot(V, dSOC_dV_to_fit, 'r');
    hold off;
    end
    
    % 7. Subtract the estimated sigmoid from the derivative data to
    % continue estimation of other peaks
    dSOC_dV_to_fit = dSOC_dV_to_fit - dXj_dU;
end

MSMR_parameters = [U0_initial(:), Q_initial(:), w_initial(:)];
MSMR_parameters = sortrows(MSMR_parameters, 1);
disp(MSMR_parameters)
num_reactions = size(MSMR_parameters, 1);

%% Optional: Generate plots for initial parameters
% Initialise arrays for Xj and dXj_dU

Xj = cell(num_reactions, 1);
dXj_dU = cell(num_reactions, 1);

for i = 1:num_reactions
    % Get parameters for this reaction
    U0 = MSMR_parameters(i, 1);
    Xj_tot = MSMR_parameters(i, 2);
    omega = MSMR_parameters(i, 3);
    
    % Get function handles for the occupancy and differential capacity
    [Xj{i}, dXj_dU{i}] = individual_reactions_function(U0, Xj_tot, omega, T);
end

X = @(x) Xj{1}(x) + Xj{2}(x) + Xj{3}(x) + Xj{4}(x);
dX_dU = @(x) dXj_dU{1}(x) + dXj_dU{2}(x) + dXj_dU{3}(x) + dXj_dU{4}(x);

if generate_plot
figure;
subplot(2, 1, 1);
plot(V, dSOC_dV, 'r');
hold on;
plot(V, dX_dU(V), 'b:.');  % Plot Voltage vs. dX/dU
xlabel('Voltage (V)');
ylabel('dSOC/dU');
legend;
hold off;

subplot(2, 1, 2);
plot(V, S, 'r');
hold on;
plot(V, X(V), 'b:.');  % Plot Voltage vs. X
xlabel('Voltage (V)');
ylabel('SOC');
legend;
hold off;
end

% Initial Guess for MSMR Parameters
initialParams = MSMR_parameters;  % Assuming MSMR_parameters is your initial guess

% Set bounds
A = [];
b = [];
Aeq = [];
beq = [];
% lb = zeros(size(initialParams)); % Lower bounds (if you have any) % Make it positive
lb = [2.5, 0, 0; 2.5, 0, 0; 2.5, 0, 0; 2.5, 0, 0];
% ub = []
ub = [4.2, Inf, Inf; 4.2, Inf, Inf; 4.2, Inf, Inf; 4.2, Inf, Inf];
nonlcon = [];

% Set options
tol = 1e-16;
options = optimoptions('fmincon','Algorithm','sqp','Display','none', ... 
                       'FunctionTolerance',tol,'StepTolerance',tol);

% Reset warning
warning('on','MATLAB:ode15s:IntegrationTolNotMet');

% Run optimisation
[MSMR_parameters, ~] = fmincon(@(p) msmrObjective(p, V, S, real_dSOC_dV), ...
                            initialParams, [], [], [], [], lb, ub, nonlcon, options);


% Display optimised parameters
disp('Optimised MSMR Parameters:');
disp(MSMR_parameters);
RMSE = msmrObjective(MSMR_parameters, V, S, real_dSOC_dV);
fprintf('RMSE: %.12f\n', RMSE);


%% Generate plots
% Create a range of voltages (U) to compute SOC
U = V;  % Voltage range (example)
num_reactions = size(MSMR_parameters, 1);

% Initialise arrays for Xj and dXj_dU
Xj = zeros(num_reactions, length(U));
dXj_dU = zeros(num_reactions, length(U));

for i = 1:num_reactions
    U0 = MSMR_parameters(i, 1);
    Xj_tot = MSMR_parameters(i, 2);
    omega = MSMR_parameters(i, 3);
    [Xj(i, :), dXj_dU(i, :)] = individual_reactions(U, U0, Xj_tot, omega, T);  % potential-dependent lithium occupancy
end

% Sum over all reactions
X = sum(Xj);
dX_dU = sum(dXj_dU);

figure;
subplot(2, 1, 1);
plot(V, real_dSOC_dV, 'r');
hold on;
% plot(V, dSOC_dV, V(TF), dSOC_dV(TF), 'g*');
plot(V, dX_dU, 'b:.');  % Plot Voltage vs. dX/dU
xlabel('Voltage (V)');
ylabel('dSOC/dU');
legend;
hold off;

subplot(2, 1, 2);
plot(V, S, 'r');
hold on;
plot(V, X, 'b:.');  % Plot Voltage vs. X
xlabel('Voltage (V)');
ylabel('SOC');
legend;
hold off;

reset_path;
ModelName = 'OCV_MSMR';
Estimator = 'PEM';
addpath(genpath(strcat('./BatEst/Code/Models/',ModelName)));
addpath(genpath(strcat('./BatEst/Code/Methods/',Estimator)));
end


function rmse = msmrObjective(MSMR_params, V, S, dSOC_dV)
    % Create a range of voltages (U) to compute SOC
    num_reactions = size(MSMR_params, 1);
    T = 25+273.15;

    % Initialise arrays for Xj and dXj_dU
    Xj = zeros(num_reactions, length(V));
    dXj_dU = zeros(num_reactions, length(V));

    for i = 1:num_reactions
        U0 = MSMR_params(i, 1);
        Xj_tot = MSMR_params(i, 2);
        omega = MSMR_params(i, 3);
        [Xj(i, :), dXj_dU(i, :)] = individual_reactions(V, U0, Xj_tot, omega, T);  % potential-dependent lithium occupancy
    end

    % Sum over all reactions
    X = sum(Xj);
    dX_dU = sum(dXj_dU);

    X_rmse = sqrt(mean((X - S).^2)) / abs(mean(S));
    dX_dU_rmse = sqrt(mean(((dX_dU - dSOC_dV)).^2)) / abs(mean(dSOC_dV));

    X_MAE = mean(abs(X - S)) / abs(mean(S));
    dX_dU_MAE = mean(abs(dX_dU - dSOC_dV)) / abs(mean(dSOC_dV));

    rmse = dX_dU_rmse + dX_dU_MAE + X_MAE;
    rmse = dX_dU_rmse;
end

function crossing = findCrossing(V, dSOC_dV, peakLoc, targetValue, direction)
    % Finds the crossing point where dSOC_dV crosses the target value.
    % V: Voltage array
    % dSOC_dV: Differential SOC array
    % peakLoc: Location of the peak
    % targetValue: The dSOC_dV value (halfMax or threequarterMax) to find crossings for
    % direction: 'left' for searching to the left, 'right' for searching to the right
    % Returns the voltage at the crossing point
    crossing = NaN; % Default value if no crossing is found
    if strcmp(direction, 'left')
        % dSOC_dV(i) > dSOC_dV(i+1)
        for i = peakLoc-1:-1:3
            if dSOC_dV(i) < dSOC_dV(i+1) && dSOC_dV(i-1) < dSOC_dV(i) && dSOC_dV(i-2) < dSOC_dV(i)
                break;
            elseif (dSOC_dV(i)-targetValue)*(dSOC_dV(i+1)-targetValue) < 0  
                % Linear interpolation to find the exact crossing point
                crossing = interp1([dSOC_dV(i), dSOC_dV(i+1)], [V(i), V(i+1)], targetValue);
                break;
            end
        end

    else % 'right'
        % dSOC_dV(i) < dSOC_dV(i+1)
        for i = peakLoc:1:length(dSOC_dV)-3
            if dSOC_dV(i) > dSOC_dV(i+1) && dSOC_dV(i+1) > dSOC_dV(i+2) && dSOC_dV(i+1) > dSOC_dV(i+3)
                break;
            elseif (dSOC_dV(i)-targetValue)*(dSOC_dV(i+1)-targetValue) < 0  
                % Linear interpolation to find the exact crossing point
                crossing = interp1([dSOC_dV(i), dSOC_dV(i+1)], [V(i), V(i+1)], targetValue);
                break;
            end
        end
    end

end

