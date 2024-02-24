function [MSMR_parameters, RMSE] = get_MSMR_params(OCP_filename)
generate_plot = true;

%% 0. Settings
% Define essential constants
T = 298.15; % temp
Rg = 8.314472;  % gas constant (J mol-1 K-1)
F = 96487;     % Faraday's constant (C mol-1)
f = F / (Rg * T);  % F/RT

% Get V, dSOC_dV data
[V, S, dSOC_dV, real_dSOC_dV] = load_dSOC_dV(OCP_filename);

%% 1. Get initial estimates
% Initialise the parameter vectors
Q_initial = zeros(1, 4);
U0_initial = zeros(1, 4);
w_initial = zeros(1, 4);

dSOC_dV_to_fit = dSOC_dV;

% Loop over the number of sigmoids to fit
for j = 1:4
    % 1) Get the first peak
    % Invert the signal to find troughs as peaks
    [pks, locs, w, p] = findpeaks(-dSOC_dV_to_fit, 'SortStr', 'descend', 'NPeaks', 5);

    peakLocs = [locs', p'.*w'.*pks', p', w', pks'];
    peakLocs = sortrows(peakLocs, 2, 'descend');
    peakLoc = peakLocs(1);

    V_peak = V(peakLoc);
    dSOC_dV_peak = dSOC_dV_to_fit(peakLoc);

    % Use the location of the first peak as U0
    U0_initial(j) = V_peak;

    % 2) Compute w_initial
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

    % Adjust position of U0_initial
    U0_initial(j) = (U0_adjustment + U0_initial(j)) * 1/2;

    % 3) compute Q_initial
    Q_initial(j) = dSOC_dV_peak * (-4) * w_initial(j) / f;
 
    % 4) Get estimated sigmoid curve
    [Xj, dXj_dU] = individual_reactions(V, U0_initial(j), Q_initial(j), w_initial(j), T); 

    if generate_plot
    figure;
    plot(V, dXj_dU, '--', 'DisplayName', 'dXj_dU');
    hold on;
    plot(V, dSOC_dV_to_fit, 'r', 'DisplayName', 'dSOC_dV_to_fit');
    legend;
    hold off;
    end
    
    % 5) Subtract the estimated sigmoid from the derivative data to
    % continue estimation of other peaks
    dSOC_dV_to_fit = dSOC_dV_to_fit - dXj_dU;
end

MSMR_parameters = [U0_initial(:), Q_initial(:), w_initial(:)];
MSMR_parameters = sortrows(MSMR_parameters, 1);
disp(MSMR_parameters)
num_reactions = size(MSMR_parameters, 1);

%% 2. Generate plots for initial parameters
[x, dx_du] = total_reaction(MSMR_parameters);

if generate_plot
figure;
% Plot Voltage vs. dX/dU
subplot(2, 1, 1);
plot(V, real_dSOC_dV, 'r', 'DisplayName', 'Real data');
hold on;
plot(V, dx_du(V), 'b:.', 'DisplayName', 'Estimated function');
xlabel('Voltage (V)');
ylabel('dSOC/dU');
legend;
hold off;

% Plot Voltage vs. X
subplot(2, 1, 2);
plot(V, S, 'r', 'DisplayName', 'Real data');
hold on;
plot(V, x(V), 'b:.', 'DisplayName', 'Estimated function');
xlabel('Voltage (V)');
ylabel('SOC');
legend;
hold off;
end

%% 3. Use Fmincon to compute accurate parameters
% Initial guess for MSMR Parameters
initialParams = MSMR_parameters;  % Assuming MSMR_parameters is your initial guess

% Set bounds
A = [];
b = [];
Aeq = [];
beq = [];
% lb = zeros(size(initialParams)); % Make it positive
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


%% 4. Generate final plots
[x, dx_du] = total_reaction(MSMR_parameters);

if generate_plot
figure;
% Plot Voltage vs. dX/dU
subplot(2, 1, 1);
plot(V, real_dSOC_dV, 'r', 'DisplayName', 'Real data');
hold on;
plot(V, dx_du(V), 'b:.', 'DisplayName', 'Estimated function');
xlabel('Voltage (V)');
ylabel('dSOC/dU');
legend;
hold off;

% Plot Voltage vs. X
subplot(2, 1, 2);
plot(V, S, 'r', 'DisplayName', 'Real data');
hold on;
plot(V, x(V), 'b:.', 'DisplayName', 'Estimated function');
xlabel('Voltage (V)');
ylabel('SOC');
legend;
hold off;
end

reset_path;
ModelName = 'OCV_MSMR';
Estimator = 'PEM';
addpath(genpath(strcat('./BatEst/Code/Models/',ModelName)));
addpath(genpath(strcat('./BatEst/Code/Methods/',Estimator)));
end


function rmse = msmrObjective(MSMR_parameters, V, S, dSOC_dV)
    [x, dx_du] = total_reaction(MSMR_parameters);
    X = x(V);
    dX_dU = dx_du(V);

    X_rmse = sqrt(mean((X - S).^2)) / abs(mean(S));
    dX_dU_rmse = sqrt(mean(((dX_dU - dSOC_dV)).^2)) / abs(mean(dSOC_dV));

    X_MAE = mean(abs(X - S)) / abs(mean(S));
    dX_dU_MAE = mean(abs(dX_dU - dSOC_dV)) / abs(mean(dSOC_dV));

    % rmse = dX_dU_rmse + X_rmse;
    rmse = dX_dU_rmse + 0.2 * X_rmse;
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

