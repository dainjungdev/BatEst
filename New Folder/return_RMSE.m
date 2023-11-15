function [RMSE_mV] = return_RMSE(Target,params,true_sol,pred_sol)
% A function to compute the root mean square error (RMSE) between the data
% and the simulation and add these values to the params structure.

RMSE_mV = sqrt(sum((true_sol.ysol(:,1)-pred_sol.ysol(:,1)).^2) ...
                        /length(pred_sol.tsol))*1e3;
% disp(RMSE_mV)
% 
% RMSE_Ts = params.Trng* ...
%                 sqrt(sum((true_sol.ysol(:,2)-pred_sol.ysol(:,2)).^2) ...
%                      /length(pred_sol.tsol)); % [C or K]
% disp(RMSE_Ts)

% verbose = params.verbose;
% 
% % Compute RMSE
% if true
%     [RMSE_mV, RMSE_Ts] = deal(0);
%     if length(true_sol.tsol)==length(pred_sol.tsol)
%         RMSE_mV = sqrt(sum((true_sol.ysol(:,1)-pred_sol.ysol(:,1)).^2) ...
%                         /length(pred_sol.tsol))*1e3; % [mV]
%         if verbose
%             disp(['The voltage RMSE is ' num2str(RMSE_mV) ' mV.']);
%         end
%         if size(true_sol.ysol,2)==2 && size(pred_sol.ysol,2)==2
%             RMSE_Ts = params.Trng* ...
%                 sqrt(sum((true_sol.ysol(:,2)-pred_sol.ysol(:,2)).^2) ...
%                      /length(pred_sol.tsol)); % [C or K]
%             if verbose
%                 disp(['The surface temperature RMSE is ' num2str(RMSE_Ts) ' K.']);
%             end
%         end
%     end
% end

end