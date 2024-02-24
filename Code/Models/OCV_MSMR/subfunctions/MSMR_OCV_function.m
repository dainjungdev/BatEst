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

[x, dx_du] = total_reaction(MSMR_parameters);
X = x(U);

% Flip X and prepend 0
X = [0, 1-X];
U = [1, U]; 

spl = pchip(X, U);
OCV = @(x, nu, miu) ppval(spl, x);

% SOC vs Voltage
if plot_model
    x = [0, logspace(-3, -0.3, 100)];
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