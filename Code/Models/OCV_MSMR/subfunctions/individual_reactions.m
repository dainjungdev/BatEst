function [xj, dxjdu] = individual_reactions(U, U0, Xj, w, T)
    % Calculates the fractional occupancy and differential capacity for an individual reaction.

    % Solves the individual reaction within an electrode. The Xj parameter can be substituted by
    % a capacity parameter (Qj) instead, and can be rescaled back to the thermodynamic factor by
    % that same capacity factor.
    % 
    % Parameters:
    % 
    % U: (float, 1-D array) Potential to be calculated
    % U0: (float) Standard electrode potential for reaction j
    % Xj: (float) Maximum fractional occupancy for reaction j (intensive) or maximum capacity of reaction j (extensive)
    % w: (float) Thermodynamic factor for reaction j
    % T: (float) temperature
    % 
    % Returns:
    % xj: (1-D array) Fractional occupancy (or capacity if extensive)
    % dxjdu: (1-D array) Differential capacity

    % Constants
    Rg = 8.314472;  % gas constant (J mol-1 K-1)
    F = 96487;     % Faraday's constant (C mol-1)
    T = 298.15;
    f = F / (Rg * T);  % F/RT

    % Fractional occupancy (xj)
    xj = Xj ./ (1 + exp(f .* (U - U0) ./ w));

    % Differential capacity (dxjdu)
    try
        dxjdu = (-Xj ./ w) .* ((f .* exp(f .* (U - U0) ./ w)) ./ (1 + exp(f .* (U - U0) ./ w)) .^ 2);
    catch
        dxjdu = zeros(size(U));  % Approximate as zero in case of Overflow Error
    end
end