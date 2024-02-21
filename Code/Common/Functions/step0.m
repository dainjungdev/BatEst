function [Model, params] = step0(ModelName,j,new_params)
% This step defines the settings, parameters and model functions.

% if strcmp(ModelName,'OCV_MSMR')
%     info = new_params.info;
%     new_params = [];
% end

% Define the parameters
params = set_parameters(j);

% If there is an input structure, overwrite with new parameters
if isstruct(new_params)
    params = convert_params(params,new_params);
end

% if strcmp(ModelName,'OCV_MSMR')
% params.info = info;
% end
% Define subfunctions
params = subfunctions(params);

% Define the model
[Model, params] = set_model(ModelName,params,j);


end
