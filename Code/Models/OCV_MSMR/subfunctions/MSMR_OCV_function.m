function OCV = MSMR_OCV_function(params)
% Unpack parameters
[OCP_filename, plot_model] = ...
    struct2array(params, {'OCP_filename','plot_model'});

[MSMR_parameters] = get_MSMR_parameters(2, 3, 1);

reset_path;
ModelName = 'OCV_MSMR';
Estimator = 'PEM';
addpath(genpath(strcat('./BatEst/Code/Models/',ModelName)));
addpath(genpath(strcat('./BatEst/Code/Methods/',Estimator)));

T = 25+273.15;   % ambient temperature (K)
Faraday = 96487;     % Faraday's constant (C mol-1)
Rg = 8.314472;       % gas constant (J mol-1 K-1)
f = Faraday / (Rg * T); % F/RT

U = linspace(2.5, 4.2, 10000);  % Voltage range (example)
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

% Create a spline representation
spl = spline(1-X, U);

OCV = @(x, nu, miu) max(2.5, min(4.2, ppval(spl, x)));

% Get MSMR_parameters(estimation from OCV data)
% Define OCV function from MSMR parameters
% Use the new OCV as OCV

% OCV_filename = OCP_filename{1};
% % Define the open-circuit voltage
% if endsWith(OCV_filename,'.parquet','IgnoreCase',true)
%     OCV = load_OCP(OCV_filename);
% elseif endsWith(OCV_filename,'.csv','IgnoreCase',true) ...
%         || endsWith(OCV_filename,'.mat','IgnoreCase',true)
%     error(['Please convert the OCP data into the Parquet file format.' ...
%            'See the DATA_PREP_GUIDE for more information.']);
% else % filename is a function
%     OCV = eval(OCV_filename);
% end

% % Define parameters for MSMR model
%     %% Single Reaction
%     % % Testing with example parameters for each reaction j
%     % U0 = 3.5;   % Standard potential for reaction j (Example)
%     % Qj_tot = 1.8;  % Total insertion capacity (Example)
%     % omega = 1;  % if Nernst equilibrium, omega = 1 (Example)
% 
%     %% Multiple Reaction
%     % Reference: Thermodynamic Model for Substitutional Materials: Application 
%     % to Lithiated Graphite, Spinel Manganese Oxide, Iron Phosphate, and Layered
%     % Nickel-Manganese-Cobalt Oxide" by Verbrugge et al., (2017)
%     % U0, Xj_tot, omega
% 
%     Faraday = params.Faraday;  % Faraday's constant
%     Rg = params.Rg;  % Gas Constant
%     T = params.Tamb; % Assuming Tamb is a temperature in Kelvin
%     f = Faraday / (Rg * T); % F/RT
% 
%     % Create a range of voltages (U) to compute SOC
%     U = linspace(3.4, 5.0, 200);  % Voltage range (example)
%     num_reactions = size(MSMR_parameters, 1);
% 
%     % Initialise arrays for Xj and dXj_dU
%     Xj = zeros(num_reactions, length(U));
%     dXj_dU = zeros(num_reactions, length(U));
% 
%     for i = 1:num_reactions
%         U0 = MSMR_parameters(i, 1);
%         Xj_tot = MSMR_parameters(i, 2);
%         omega = MSMR_parameters(i, 3);
%         [Xj(i, :), dXj_dU(i, :)] = individual_reactions(U, U0, Xj_tot, omega, T);  % potential-dependent lithium occupancy
%     end
% 
%     % Sum over all reactions
%     X = sum(Xj);
%     dX_dU = sum(dXj_dU);
% 
%     try
%     % Create a spline representation
%     spl = spline(X, U);
% 
%     % Define the OCV function using ppval
%     % Define OCV as a subfunction
%     OCV_o = @(x) ppval(spl, x);
%     OCV = @(SOC,nu,miu) OCV_o(miu-nu*SOC);
%     catch
%         plot_model_2 = false;
%     end
% 
% %% need function V = OCV(x, c1, c2, c3)
% 
%    % We need unknown parameters instead of miu and nu
% 
%     % Plot OCV if plot_model is true
%     % Voltage vs. Xj, Voltage vs dXj/dU
%     if plot_model
%         subplot(2, 2, 1);
%         hold on;
%         plot(U, Xj, '.');  % Plot Voltage vs. Xj
%         xlabel('Voltage (V)');
%         ylabel('Xj');
%         legend;
%         hold off;
% 
%         subplot(2, 2, 2);
%         hold on;
%         plot(U, dXj_dU, '.');  % Plot Voltage vs. dXj/dU
%         xlabel('Voltage (V)');
%         ylabel('dXj/dU');
%         legend;
%         hold off;
% 
%         subplot(2, 2, 3);
%         hold on;
%         plot(U, X, 'b:.');  % Plot Voltage vs. X
%         xlabel('Voltage (V)');
%         ylabel('X');
%         hold off;
% 
%         subplot(2, 2, 4);
%         hold on;
%         plot(U, dX_dU, 'b:.');  % Plot Voltage vs. dX/dU
%         xlabel('Voltage (V)');
%         ylabel('dX/dU');
%         hold off;
%     end


% SOC vs Voltage
if plot_model
    x = [0, logspace(-3,-0.3,100)];
    x = unique([x,1-x]);
    figure; hold on;
    plot(x,OCV(x),'b:.');
    xlabel('State of charge (SOC)');
    ylabel('Voltage (V)')
end
end