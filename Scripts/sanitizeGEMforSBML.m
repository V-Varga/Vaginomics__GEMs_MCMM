function model = sanitizeGEMforSBML(model)
% sanitizeGEMforSBML  Clean a RAVEN draft GEM for SBML export
%
%   model = sanitizeGEMforSBML(model)
%
%   Purpose:
%       Ensures that a draft genome-scale metabolic model (GEM) generated 
%       with RAVEN can be exported to SBML in a way that is compatible with 
%       COBRApy, GLPK, and other SBML-compliant software. 
%
%   Specifically:
%       • Metabolite, reaction, and compartment IDs are restricted to 
%         alphanumeric characters and underscores (_).
%       • Names for metabolites, reactions, and compartments are also 
%         restricted to alphanumeric + underscore.
%       • Chemical formulas containing parentheses, the character "R" 
%         (indicating generic groups), or polymer repeats are removed, 
%         since these are not SBML-compliant.
%
%   Input:
%       model - RAVEN model struct (draft GEM)
%
%   Output:
%       model - sanitized model struct safe for SBML export
%
%   See also: prepareModelForExport, exportModel

    % --- Helper function to clean strings ---
    % Replace any character not [a-zA-Z0-9_] with an underscore
    sanitizeStr = @(s) regexprep(s, '[^a-zA-Z0-9_]', '_');

    %--------------------------------------------------------------
    % Sanitize metabolite IDs, names, and formulas
    %--------------------------------------------------------------
    for i = 1:length(model.mets)
        % Sanitize metabolite IDs (unique identifiers)
        model.mets{i} = sanitizeStr(model.mets{i});

        % Sanitize metabolite names (if present)
        if isfield(model, 'metNames') && i <= length(model.metNames)
            model.metNames{i} = sanitizeStr(model.metNames{i});
        end

        % Remove problematic chemical formulas (if present)
        if isfield(model, 'metFormulas') && i <= length(model.metFormulas)
            f = model.metFormulas{i};
            % If formula contains parentheses, "R", or other non-SBML-safe
            % characters, clear the entry entirely
            if ~isempty(f) && any(contains(f, {'(', ')', 'R'}))
                model.metFormulas{i} = ''; 
            end
        end
    end

    %--------------------------------------------------------------
    % Sanitize reaction IDs and names
    %--------------------------------------------------------------
    for i = 1:length(model.rxns)
        % Sanitize reaction IDs
        model.rxns{i} = sanitizeStr(model.rxns{i});

        % Sanitize reaction names (if present)
        if isfield(model, 'rxnNames') && i <= length(model.rxnNames)
            model.rxnNames{i} = sanitizeStr(model.rxnNames{i});
        end
    end

    %--------------------------------------------------------------
    % Sanitize compartment IDs and names
    %--------------------------------------------------------------
    if isfield(model, 'comps')
        for i = 1:length(model.comps)
            % Sanitize compartment IDs
            model.comps{i} = sanitizeStr(model.comps{i});

            % Sanitize compartment names (if present)
            if isfield(model, 'compNames') && i <= length(model.compNames)
                model.compNames{i} = sanitizeStr(model.compNames{i});
            end
        end
    end
end
