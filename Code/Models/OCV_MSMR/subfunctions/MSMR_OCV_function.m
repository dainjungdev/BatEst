function params = MSMR_OCV_function(params)
% Unpack parameters
[OCP_filename, plot_model] = ...
    struct2array(params, {'OCP_filename','plot_model'});

MSMR_parameters = get_MSMR_parameters(OCP_filename);

T = params.Tref;  % ambient temperature (K)
Faraday = params.Faraday;     % Faraday's constant (C mol-1)
Rg = params.Rg;       % gas constant (J mol-1 K-1)
f = Faraday / (Rg * T); % F/RT

data = parquetread(OCP_filename{1});
U = linspace(min(data.Voltage_V), max(data.Voltage_V), 10000);  % Voltage range (example)

num_reactions = size(MSMR_parameters, 1);

% Initialise arrays for Xj and dXj_dU
Xj = zeros(num_reactions, length(U));
dXj_dU = zeros(num_reactions, length(U));

% Create function handles for each reaction
reactionHandles = cell(num_reactions, 1);
for i = 1:num_reactions
    % Get parameters for this reaction
    U0 = MSMR_parameters(i, 1);
    Xj_tot = MSMR_parameters(i, 2);
    omega = MSMR_parameters(i, 3);
    
    % Get function handles for the occupancy and differential capacity
    reactionHandles{i} = individual_reactions_function(U0, Xj_tot, omega, T);
end

% Evaluate the function handles over the voltage range U
for i = 1:num_reactions
    Xj(i, :) = reactionHandles{i}.xj(U);
    dXj_dU(i, :) = reactionHandles{i}.dxjdu(U);
end

% Sum over all reactions
X = sum(Xj, 1);  % Sum across rows for each voltage
dX_dU = sum(dXj_dU, 1);  % Sum across rows for each voltage

X = [0, 1-X];  % Prepend 0 to 1-X to include the origin point in the x-values
U = [0, U]; 

spl = pchip(X, U);
OCV = @(x, nu, miu) ppval(spl, x);
% % Create a spline representation
% spl = spline(1-X, U);
% OCV = @(x, nu, miu) max(0, min(4.2, ppval(spl, x))); % Ensure OCV is within the voltage range [2.5, 4.2]

% Assuming 1-X and U are your data vectors
% OCV = @(x, nu, miu) max(1.5, min(4.2, interp1(1-X, U, x, 'linear', 'extrap')));

% SOC vs Voltage
if plot_model
    x = [0, logspace(-3,-0.3,100)];
    x = unique([x,1-x]);
    figure;
    plot(x,OCV(x),'b');
    xlabel('State of charge (SOC)');
    ylabel('Voltage (V)')
    title('x vs OCV(x)');
end

params.OCV = OCV;
params.MSMR_parameters = reshape(MSMR_parameters, 1, []);

end