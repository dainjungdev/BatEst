function observe_changes
    % Read the parquet file
    myTable = parquetread('2_cycling_analysis/out/EHM_24.parquet');

    % List of varying parameters
    varying_params = {'Q', 'tau', 'b', 'Ip', 'In', 'nu', 'miu', 'Rf'};

    figure;

    % Loop through varying_params
    for i = 1:length(varying_params)
        paramName = varying_params{i};

        % Extract data for plotting
        x_data = myTable.RowN;
        y_data = myTable.(paramName); % Dynamic field reference

        % Create a plot
        plot(x_data, y_data, '-o');
        title(['Plot for ' paramName]);
        xlabel('RowN');
        ylabel(paramName);

        % Create a subplot
        subplot(4, 2, i);
        plot(x_data, y_data, '-o');
        title(['Plot for ' paramName]);
        xlabel('RowN');
        ylabel(paramName);
    end

    % Adjust layout
    set(gcf, 'Position', get(0, 'Screensize'));
end
