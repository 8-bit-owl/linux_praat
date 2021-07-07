function parameters = checkParameters(parameters,defaults)
% This function checks to see if all fields in defaults are in parameters.  Any fields not in
% parameters are copied over from defaults.
% Input:
%    parameters -- struct of parameters used in an experiment (if empty, output is defaults)
%    defaults -- struct of default parameter values.  If a field in defaults is not in parameters,
%       the field is created in parameters and initialized with the field value in defaults.
% Output:
%    parameters -- same as input, with added fields and values from defaults
% Note: this function does NOT check that parameter values are VALID.

% Check for empty parameters:
if isempty(parameters),
   parameters = defaults;
   return;
end;

% Get field names:
parameterFields = fieldnames(parameters);
defaultFields = fieldnames(defaults);

% Check to see that each default field is in parameters, add if missing:
for p=1:length(defaultFields),
   F = defaultFields{p}; % string of default field name
   if isempty(strmatch(F,parameterFields,'exact')),
      % Add field to parameters, set to default value:
      parameters = setfield(parameters,F,getfield(defaults,F));
   end;
end;

return;

% Bye!