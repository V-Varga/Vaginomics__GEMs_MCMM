function modelOut = padMissingFields(modelIn, templateModel)
% padMissingFields  Add missing fields to a GEM using a template model
%
%   modelOut = padMissingFields(modelIn, templateModel)
%
%   PURPOSE:
%       Ensures that all fields present in a template model exist in the
%       input GEM. If a field is missing, it is added with appropriate
%       default values (empty cells, zeros, or empty structs).
%
%   INPUTS:
%       modelIn        - GEM struct (RAVEN/COBRA style)
%       templateModel  - Template GEM struct defining expected fields
%
%   OUTPUT:
%       modelOut       - Updated GEM struct with missing fields padded
%
%   NOTES:
%       - Reaction-level fields (e.g., 'rxnReferences', 'rxnConfidenceScores')
%         are padded to the number of reactions.
%       - Metabolite-level fields (e.g., 'metCharges', 'metFrom') are padded
%         to the number of metabolites.
%       - Gene-level fields (e.g., 'geneFrom', 'geneMiriams') are padded
%         to the number of genes.

    % Initialize output as input
    modelOut = modelIn;

    % Loop through all fields in the template model
    for f = fieldnames(templateModel)'
        fname = f{1};
        
        % Only act if field is missing
        if ~isfield(modelOut, fname)
            switch fname
                % Reaction-level fields
                case {'rxnFrom','rules','equations','rxnReferences'}
                    modelOut.(fname) = repmat({''}, numel(modelOut.rxns), 1);
                case 'rxnConfidenceScores'
                    modelOut.(fname) = zeros(numel(modelOut.rxns),1);
                case 'grRules'
                    modelOut.(fname) = repmat({''}, numel(modelOut.rxns), 1);

                % Metabolite-level fields
                case {'metFrom'}
                    modelOut.(fname) = repmat({''}, numel(modelOut.mets), 1);
                case 'metCharges'
                    modelOut.(fname) = zeros(numel(modelOut.mets),1);

                % Gene-level fields
                case {'geneFrom','geneNames','geneShortNames'}
                    modelOut.(fname) = repmat({''}, numel(modelOut.genes), 1);
                case 'geneMiriams'
                    modelOut.(fname) = repmat({struct()}, numel(modelOut.genes), 1);

                % Default: empty cell array
                otherwise
                    modelOut.(fname) = {};
            end
            fprintf('Added %s to KEGG\n', fname);
        end
    end
end
