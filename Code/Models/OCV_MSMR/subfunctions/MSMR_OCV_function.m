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
U = linspace(min(data.Voltage_V), max(data.Voltage_V), 1000);  % Voltage range (example)
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
OCV = @(x, nu, miu) max(2.5, ppval(spl, x));

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
end

params.OCV = OCV;
params.MSMR_parameters = reshape(MSMR_parameters, 1, []);

end