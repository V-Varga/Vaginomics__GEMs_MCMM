function model = sanitizeRxnMiriams(model)
% sanitizeRxnMiriams  Ensure the rxnMiriams field is present and well-formed
%
%   model = sanitizeRxnMiriams(model)
%
%   Purpose:
%       Ensures that every reaction in the model has a properly initialized
%       MIRIAM annotation structure. This is important for SBML compliance 
%       and for maintaining consistent annotation data across reactions.
%
%   Behavior:
%       • Creates the rxnMiriams field if missing.
%       • For each reaction:
%           - If empty, assigns an empty struct with fields 'name' and 'value'.
%           - If the entry exists but is not a struct, replaces it with an empty struct.
%           - If the struct exists but lacks 'name' or 'value', adds those fields.
%
%   Input:
%       model - RAVEN/COBRA model struct
%
%   Output:
%       model - updated model with sanitized rxnMiriams field

    % --- Ensure the rxnMiriams field exists ---
    if ~isfield(model,'rxnMiriams')
        model.rxnMiriams = cell(numel(model.rxns),1);
    end

    % --- Iterate over all reactions ---
    for i = 1:numel(model.rxns)
        if isempty(model.rxnMiriams{i})
            % Entry completely missing → create empty struct
            model.rxnMiriams{i} = struct('name',{{}},'value',{{}});
        elseif ~isstruct(model.rxnMiriams{i})
            % Wrong type → replace with empty struct
            model.rxnMiriams{i} = struct('name',{{}},'value',{{}});
        else
            % Struct exists → ensure both required fields exist
            if ~isfield(model.rxnMiriams{i},'name')
                model.rxnMiriams{i}.name = {};
            end
            if ~isfield(model.rxnMiriams{i},'value')
                model.rxnMiriams{i}.value = {};
            end
        end
    end

    % --- Status message ---
    fprintf('Sanitized rxnMiriams for %d reactions\n', numel(model.rxns));
end
