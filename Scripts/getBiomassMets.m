function biomassMets = getBiomassMets(model, biomassRxnID)
% getBiomassMets  Extract metabolite IDs involved in the biomass reaction
%
%   biomassMets = getBiomassMets(model, biomassRxnID)
%
%   PURPOSE:
%       Returns all metabolites (KEGG IDs) that participate in the biomass reaction
%       of a GEM (Genome-scale Metabolic Model) in RAVEN/COBRA format.
%
%   INPUTS:
%       model        - struct, GEM in RAVEN or COBRA style
%       biomassRxnID - string, optional reaction ID for biomass. 
%                      If omitted, the reaction with objective coefficient c==1
%                      will be used.
%
%   OUTPUT:
%       biomassMets  - cell array of metabolite IDs involved in biomass
%
%   NOTES:
%       Displays the reaction ID and number of metabolites for convenience.

    % --- Step 1: Identify the biomass reaction ---
    if nargin < 2 || isempty(biomassRxnID)
        rxnIdx = find(model.c ~= 0); % use objective function coefficient
        if isempty(rxnIdx)
            error('No biomass reaction found (c vector is all zero).');
        elseif numel(rxnIdx) > 1
            warning('Multiple biomass reactions found, using the first one.');
            rxnIdx = rxnIdx(1);
        end
    else
        rxnIdx = find(strcmp(model.rxns, biomassRxnID));
        if isempty(rxnIdx)
            error('Biomass reaction %s not found in model.', biomassRxnID);
        end
    end

    % --- Step 2: Extract the stoichiometry column for the reaction ---
    stoich = model.S(:, rxnIdx);

    % --- Step 3: Identify participating metabolites (non-zero stoichiometry) ---
    metIdx = find(stoich ~= 0);

    % --- Step 4: Return metabolite IDs ---
    biomassMets = model.mets(metIdx);

    % --- Display reaction information ---
    fprintf('Biomass reaction: %s\n', model.rxns{rxnIdx});
    fprintf('Number of metabolites in biomass reaction: %d\n', numel(biomassMets));
    disp(biomassMets);
end
