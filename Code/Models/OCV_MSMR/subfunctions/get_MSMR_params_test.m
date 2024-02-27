function [MSMR_parameters, RMSE] = get_MSMR_params_test(OCP_filename)
generate_plot = true;

%% 0. Settings
% Define essential constants
T = 298.15; % temp
Rg = 8.314472;  % gas constant (J mol-1 K-1)
F = 96487;     % Faraday's constant (C mol-1)
f = F / (Rg * T);  % F/RT

% Get V, dSOC_dV data
[V, S, dSOC_dV, real_dSOC_dV, V_ext, S_ext, dSOC_dV_ext] = load_dSOC_dV_edit(OCP_filename);
% V = [V_ext, V];
% S = [ones(1, length(V_ext)), S];
% dSOC_dV = [dSOC_dV_ext, dSOC_dV];
% real_dSOC_dV = [dSOC_dV_ext, real_dSOC_dV];


V = [V_ext V];
S = [S_ext S];
dSOC_dV = [dSOC_dV_ext dSOC_dV];
real_dSOC_dV = [dSOC_dV_ext real_dSOC_dV];


%% 1. Get initial estimates
% Initialise the parameter vectors
Q_initial = zeros(1, 2);
U0_initial = zeros(1, 2);
w_initial = zeros(1, 2);

dSOC_dV_to_fit = dSOC_dV;

% Define peak location bounds
peakLoc_bounds = [800 1500; 1200 1600; 1400 1800; 1700 2000];
peakLoc_bounds = peakLoc_bounds + length(V_ext);

% Reorder bounds into desired sequence
sequence = [3, 4]; % Define the sequence of rows you want
peakLoc_limits = peakLoc_bounds(sequence, :);

% Loop over the number of sigmoids to fit
for j = 1:2
    % 1) Get the first peak
    % Invert the signal to find troughs as peaks
    [pks, locs, w, p] = findpeaks(-dSOC_dV_to_fit, 'SortStr', 'descend', 'NPeaks', 7);

    % Combine the peak information into a single matrix for easier processing
    peakInfo = [locs', (p'.^0.9).*(w'.*pks'), p', w', pks'];
    peakInfo = sortrows(peakInfo, 2, 'descend');
    
    % Filter the peaks based on their locations (1200 to 1800)
    peakLoc_limit = peakLoc_limits(j,:);
    filteredPeakInfo = peakInfo(peakInfo(:,1) >= peakLoc_limit(1) & peakInfo(:,1) <= peakLoc_limit(2), :);
    
    if ~isempty(filteredPeakInfo)
        peakLoc = filteredPeakInfo(1);
    else
        peakLoc = peakInfo(1);
    end

    V_peak = V(peakLoc);
    dSOC_dV_peak = dSOC_dV_to_fit(peakLoc);

    % Use the location of the first peak as U0
    U0_initial(j) = V_peak;

    % 2) Compute w_initial
    % range = [linspace(0.9, 0.5, 100), linspace(0.5, 0.05, 200)];
    range = [linspace(0.99, 0.01, 200)] .^ 1.2;
    % range = 0.995:-0.001:0.005;

    U0_adjustment = U0_initial(j);
    counter = 0;
    for i = 1:length(range)
        k = range(i);
        ref_dSOC_dV = dSOC_dV_peak * k;
        V_ref_left = findCrossing(V, dSOC_dV_to_fit, peakLoc, ref_dSOC_dV, 'left');
        V_ref_right = findCrossing(V, dSOC_dV_to_fit, peakLoc, ref_dSOC_dV, 'right');

        if ~isnan(V_ref_left) && ~isnan(V_ref_right)
            ref_distance_left = V_peak - V_ref_left;
            ref_distance_right = V_ref_right - V_peak;
            ref_distance = (V_ref_right - V_ref_left) * 1/2;
            % ref_distance = 2/3 * min(ref_distance_left, ref_distance_right) + 1/3 * max(ref_distance_left, ref_distance_right);
            w = ref_distance * f / log(2/k - 1 + 2*sqrt(1/(k*k) - 1/k));
            w_initial(j) = (w_initial(j) * (i - 1) + w) / i;
    
            ref_halfpoint = (V_ref_right + V_ref_left) * 1/2;
            % if ref_distance_left < ref_distance_right
            %     ref_halfpoint = V_ref_right * 9/20 + V_ref_left * 11/20;
            % else
            %     ref_halfpoint = V_ref_right * 11/20 + V_ref_left * 9/20;
            % end
            
            % if i < 0.9 * length(range)
            % U0_adjustment = (U0_adjustment * (i - 1) + ref_halfpoint) / i;
            % end

        elseif ~isnan(V_ref_left)
            ref_distance = V_peak - V_ref_left;
            w = ref_distance * f / log(2/k - 1 + 2*sqrt(1/(k*k) - 1/k));
            counter = counter + 1;
            w_initial(j) = (w_initial(j) * (i - 1/2 - 1/2 * counter) + w * 1/2) / (i - 1/2 * counter);
            if counter >= min(length(range) - i, max(1,i*1/3)) break; end

        elseif ~isnan(V_ref_right)
            ref_distance = V_ref_right - V_peak;
            w = ref_distance * f / log(2/k - 1 + 2*sqrt(1/(k*k) - 1/k));
            counter = counter + 1;
            w_initial(j) = (w_initial(j) * (i - 1/2 - 1/2 * counter) + w * 1/2) / (i - 1/2 * counter);
            if counter >= min(length(range) - i, max(1,i*1/3)) break; end

        else
            break;
        end

    end

    % Adjust w_initial for error
    w_initial(j) = w_initial(j) ;

    % Adjust position of U0_initial
    U0_initial(j) = (U0_adjustment*(i) + U0_initial(j)*length(range)/2) / (i+length(range)/2);

    % 3) compute Q_initial
    Q_initial(j) = dSOC_dV_peak * (-4) * w_initial(j)/ f; 

 
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

old_MSMR_parameters = [U0_initial(:), Q_initial(:), w_initial(:)];
% MSMR_parameters = sortrows(MSMR_parameters, 1);
disp(old_MSMR_parameters)


%% 2. Generate plots for initial parameters
[x, dx_du] = total_reaction(old_MSMR_parameters);

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








% Invert the signal to find troughs as peaks
[pks, locs, w, p] = findpeaks(-dSOC_dV_to_fit, 'SortStr', 'descend', 'NPeaks', 5);

% Combine the peak information into a single matrix for easier processing
peakInfo = [locs', (p'.^0.9).*(w'.*pks'), p', w', pks'];
peakInfo = sortrows(peakInfo, 2, 'descend');

peakLoc = peakInfo(1);
V_peak = V(peakLoc);
dSOC_dV_peak = dSOC_dV_to_fit(peakLoc);


Q_left = 1 - sum(Q_initial);
w_left = Q_left * f / (-4 * dSOC_dV_peak);
MSMR_parameters = [V_peak,Q_left*1/3,w_left; V_peak, Q_left*2/3, w_left];






%% 3. Use Fmincon to compute accurate parameters
% Initial guess for MSMR Parameters
initialParams = MSMR_parameters;  % Assuming MSMR_parameters is your initial guess

% Set bounds
A = [];
b = [];
A = [0, 0, 1, 1, 0, 0];
b = Q_left;
Aeq = []; 
beq = []; 
lb = [];
ub = Inf(size(initialParams));
nonlcon = [];

% Set options
tol = 1e-16;
options = optimoptions('fmincon','Algorithm','sqp','Display','none', ... 
                       'FunctionTolerance',tol,'StepTolerance',tol);

% Reset warning
warning('on','MATLAB:ode15s:IntegrationTolNotMet');

% Run optimisation
[MSMR_parameters, ~] = fmincon(@(p) msmrObjective(p, V, S, dSOC_dV_to_fit, 'dxdu_rmse'), ...
                            initialParams, A, b, Aeq, beq, lb, ub, nonlcon, options);

MSMR_parameters = [old_MSMR_parameters ; MSMR_parameters];

% Display optimised parameters
disp('Optimised MSMR Parameters:');
disp(MSMR_parameters);
RMSE = msmrObjective(MSMR_parameters, V, S, real_dSOC_dV, 'dxdu_rmse');
fprintf('RMSE: %.12f\n', RMSE);







%% 3. Use Fmincon to compute accurate parameters
% Initial guess for MSMR Parameters
initialParams = MSMR_parameters;  % Assuming MSMR_parameters is your initial guess

% Set bounds
A = [];
b = [];
A = [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0];
b = 1;
Aeq = []; % Aeq = [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0];
beq = []; % beq = 1;
% lb = zeros(size(initialParams)); % Make it positive
lb = [];
% lb(:,1) = initialParams(:,1) * 0.95;
% lb = [2.5, 0, 0; 2.5, 0, 0; 2.5, 0, 0; 2.5, 0, 0];
% lb = [3.1, 0, 0; 3.4, 0, 0; 3.7, 0, 0; 3.95, 0, 0];
% ub = [];
ub = Inf(size(initialParams));
% ub(:,1) = initialParams(:,1) * 1.05;
% ub(:,1) = initialParams(:,1) * 1.2;
% ub = [4.2, Inf, Inf; 4.2, Inf, Inf; 4.2, Inf, Inf; 4.2, Inf, Inf];
% ub = [3.5, Inf, Inf; 3.8, Inf, Inf; 4.05, Inf, Inf; 4.2, Inf, Inf];
nonlcon = [];

% Set options
tol = 1e-16;
options = optimoptions('fmincon','Algorithm','sqp','Display','none', ... 
                       'FunctionTolerance',tol,'StepTolerance',tol);

% Reset warning
warning('on','MATLAB:ode15s:IntegrationTolNotMet');

% Run optimisation
[MSMR_parameters, ~] = fmincon(@(p) msmrObjective(p, V, S, real_dSOC_dV, 'dxdu_rmse'), ...
                            initialParams, A, b, Aeq, beq, lb, ub, nonlcon, options);


% Display optimised parameters
disp('Optimised MSMR Parameters:');
disp(MSMR_parameters);
RMSE = msmrObjective(MSMR_parameters, V, S, real_dSOC_dV, 'dxdu_rmse');
fprintf('RMSE: %.12f\n', RMSE);






%% 4. Generate final plots
[x, dx_du] = total_reaction(MSMR_parameters);

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


% reset_path;
% ModelName = 'OCV_MSMR';
% Estimator = 'PEM';
% addpath(genpath(strcat('./BatEst/Code/Models/',ModelName)));
% addpath(genpath(strcat('./BatEst/Code/Methods/',Estimator)));
end


function rmse = msmrObjective(MSMR_parameters, V, S, dSOC_dV, metric)
    [x, dx_du] = total_reaction(MSMR_parameters);
    X = x(V);
    dX_dU = dx_du(V);

    % X_rmse = sqrt(mean((X - S).^2)) / abs(mean(S));
    % dX_dU_rmse = sqrt(mean(((dX_dU - dSOC_dV)).^2)) / abs(mean(dSOC_dV));

    X_rmse = sqrt(mean((X - S).^2)) / abs(mean(S));
    dX_dU_rmse = sqrt(mean(((dX_dU - dSOC_dV)).^2)) / abs(mean(dSOC_dV));
    dX_dU_cubic = nthroot(mean(abs(dX_dU - dSOC_dV).^3), 3) / abs(mean(dSOC_dV));

    X_MAE = mean(abs(X - S)) / abs(mean(S));
    dX_dU_MAE = mean(abs(dX_dU - dSOC_dV)) / abs(mean(dSOC_dV));

    % rmse = dX_dU_rmse + X_rmse;
    rmse = dX_dU_rmse;

    if strcmp(metric, ('dxdu_rmse'))
        rmse = dX_dU_rmse;
    elseif strcmp(metric, ('dxdu_MAE'))
        rmse = dX_dU_MAE;
    elseif strcmp(metric, ('dxdu_cubic'))
        rmse = dX_dU_cubic;
    end
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
        for i = peakLoc-2:-1:2
            if all(diff(dSOC_dV(i-1:i+1)) > 5e-05) && max(diff(dSOC_dV(i-1:i+1))) > 1e-04
                break;
            elseif (dSOC_dV(i)-targetValue)*(dSOC_dV(i+1)-targetValue) < 0  
                % Linear interpolation to find the exact crossing point
                crossing = interp1([dSOC_dV(i), dSOC_dV(i+1)], [V(i), V(i+1)], targetValue);
                break;
            end
        end

    else % 'right'
        % dSOC_dV(i) < dSOC_dV(i+1)
        for i = peakLoc+2:1:length(dSOC_dV)-1%
            if all(diff(dSOC_dV(i-1:i+1)) < -5e-05) && min(diff(dSOC_dV(i-1:i+1))) < -1e-04
                break;
            elseif (dSOC_dV(i)-targetValue)*(dSOC_dV(i+1)-targetValue) < 0  
                % Linear interpolation to find the exact crossing point
                crossing = interp1([dSOC_dV(i), dSOC_dV(i+1)], [V(i), V(i+1)], targetValue);
                break;
            end
        end
    end

end

