function reactionHandle = individual_reactions_function(U0, Xj, w, T)
    % Constants
    Rg = 8.314472;  % gas constant (J mol-1 K-1)
    F = 96487;      % Faraday's constant (C mol-1)
    
    f = F / (Rg * T);  % F/RT
    
    % Return a function handle that computes the fractional occupancy and differential capacity
    reactionHandle.xj = @(U) Xj ./ (1 + exp(f .* (U - U0) ./ w));
    reactionHandle.dxjdu = @(U) (-Xj ./ w) .* ((f .* exp(f .* (U - U0) ./ w)) ./ (1 + exp(f .* (U - U0) ./ w)) .^ 2);
end
