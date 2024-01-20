function OCV = OCV_function(params)
% This function defines the open-circuit voltage (V) as a continuous
% function of state of charge.

% Unpack parameters
[OCP_filename, plot_model] = ...
    struct2array(params, {'OCP_filename','plot_model'});

% Check if MSMR model is to be used
if isfield(params, 'MSMR') && params.MSMR
    % Define parameters for MSMR model
    % Testing with example parameters for each reaction j
    % U0 = params.U0;  % Standard potential for reaction j
    U0 = 3.5;   % Standard potential for reaction j (Example)
    % Qj_tot = params.Qj_tot;  % Total insertion capacity for reaction j
    Qj_tot = 1.8;  % Total insertion capacity (Example)
    % omega = params.omega;  % Non-ideality factor for reaction j
    omega = 1;  % (if Nernst equilibrium, omega = 1)
    Faraday = params.Faraday;  % Faraday's constant
    Rg = params.Rg;  % Gas Constant
    T = params.Tamb; % Assuming Tamb is a temperature in Kelvin
    f = Faraday / (Rg * T); % F/RT
    
    % % MSMR model for each reaction j
    % % U is the independent variable representing the voltage
    % Qj = @(U) Qj_tot ./ (1 + exp(f .* (U - U0) ./ omega));
    % 
    % % % Differential capacity for each reaction j
    % % dQj_dU = @(U) -Qj_tot .* f .* exp(f .* (U - U0) ./ omega) ./ ...
    % %              (omega .* (1 + exp(f .* (U - U0) ./ omega)).^2);
    % 
    % % % Total OCV and differential capacity
    % % OCV = @(U) sum(arrayfun(Qj, U));  % Sum of Qj(U) over all reactions
    % % dQ_dU = @(U) sum(arrayfun(dQj_dU, U));  % Sum of dQj_dU(U) over all reactions
    % 
    % % ... [rest of the MSMR implementation] ...
    % % Q is SOC, U is Voltage
    % 
    % % Plot OCV if plot_model is true
    % if plot_model
    %     % Plot open-circuit voltage based on MSMR model
    %     x = [0, logspace(-3, -0.3, 100)];  % State of charge
    %     x = unique([x, 1 - x]);
    %     y = arrayfun(Qj, x);  % Voltage for each SOC value
    %     figure; hold on;
    %     plot(x, y, 'b:+');  % Plot SOC on x-axis and Voltage on y-axis
    %     xlabel('State of charge');
    %     ylabel('Voltage (V)');
    % end

    % Create a range of voltages (U) to compute SOC
    x = linspace(2.5, 4.2, 100);  % Voltage range (example)
    Qj = individual_reactions(x, U0, Qj_tot, omega, T);  % potential-dependent lithium occupancy
    y = Qj / Qj_tot; % Calculate SOC for each voltage



    % Create a spline representation
    spl = spline(x, y);

    % Define the OCV function using ppval
    OCP = @(x) ppval(spl, x);
    OCV = @(SOC,nu,miu) OCP(miu-nu*SOC);

    % Plot OCV if plot_model is true
    if params.plot_model
        figure; hold on;
        plot(x, y, 'b:+');  % Plot Voltage vs. SOC
        xlabel('Voltage (V)');
        ylabel('State of charge');
    end


else
    if length(OCP_filename)==1
        OCV_filename = OCP_filename{1};
    
        % Define the open-circuit voltage
        if endsWith(OCV_filename,'.parquet','IgnoreCase',true)
            OCV = load_OCP(OCV_filename);
        elseif endsWith(OCV_filename,'.csv','IgnoreCase',true) ...
                || endsWith(OCV_filename,'.mat','IgnoreCase',true)
            error(['Please convert the OCP data into the Parquet file format.' ...
                   'See the DATA_PREP_GUIDE for more information.']);
        else % filename is a function
            OCV = eval(OCV_filename);
        end
        
        if plot_model
            % Plot open-circuit voltage
            x = [0, logspace(-3,-0.3,100)];
            x = unique([x,1-x]);
            figure; hold on;
            plot(x,OCV(x),'b:+');
            xlabel('State of charge');
            ylabel('Voltage (V)')
        end
        
    elseif length(OCP_filename)==2
        
        % Define the OCV as the difference between the electrode OCPs
        % UpFun = cathode potential, Unfun = anode potential
        [UnFun, UpFun] = electrode_potentials(params);
        OCV =  @(soc,nu,miu) UpFun(soc,nu,miu)-UnFun(soc);
        
    else
        error('Unexpectd OCP_filename.')
    end
end
end
