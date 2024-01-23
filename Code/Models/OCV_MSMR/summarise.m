function summary_table = summarise(out)
    % Create a new table with the selected columns
    summary_table = table(out.RowN, out.Qn ./ out.hr, out.nu, out.miu, out.S0, out.RMSE_mV, ...
                          'VariableNames', {'RowN', 'Qn/hr', 'nu', 'miu', 'S0', 'RMSE_mV'});

    % % Sort the table
    % if width(out) > 2 % Check if 'out' has more than 2 columns (e.g., 'Cell_Number' and 'Test_Number')
    %     summary_table = movevars(summary_table, 'Cell_Number', 'Before', 1);
    %     summary_table = movevars(summary_table, 'Test_Number', 'After', 'Cell_Number');
    %     summary_table = movevars(summary_table, 'RowN', 'After', 'Test_Number');
    %     summary_table = sortrows(summary_table, {'Cell_Number', 'Test_Number', 'RowN'});
    % else
    %     summary_table = movevars(summary_table, 'RowN', 'Before', 1);
    %     summary_table = sortrows(summary_table, {'RowN'});
    % end

    % % Display summary information
    % disp(['Qn = ' num2str(out.Qn(1) / hr) '*hr; % negative electrode capacity (As)']);
    % disp(['nu = ' num2str(out.nu(1)) ';       % negative-positive electrode capacity ratio (non-dim.)']);
    % disp(['miu = ' num2str(out.miu(1)) ';     % cyclable lithium-positive electrode capacity ratio (non-dim.)']);
    % disp(['S0 = ' num2str(out.S0(1)) '; % initial SOC']);
    % disp(['RMSE_mV = ' num2str(out.RMSE_mV(1))]);
    
    % Optionally, save the summary_table to a file
    % writetable(summary_table, 'summary_table.csv');
end
