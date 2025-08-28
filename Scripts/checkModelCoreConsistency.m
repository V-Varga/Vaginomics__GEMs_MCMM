function checkModelCoreConsistency(model)
% checkModelCoreConsistency  Perform basic structural checks on a COBRA/RAVEN model.
%
%   checkModelCoreConsistency(model)
%
%   PURPOSE:
%       This function performs quick structural checks on a metabolic model
%       to help debug issues (e.g., failures when exporting the model).
%       It verifies consistency between reaction, metabolite, and gene counts
%       and their corresponding matrices.
%
%   INPUTS:
%       model   - struct, COBRA/RAVEN model structure
%
%   OUTPUTS:
%       None (results are printed to the MATLAB command window)
%
%   CHECKS PERFORMED:
%       - Number of reactions, metabolites, and genes
%       - Dimensions of the stoichiometric matrix (S)
%       - Dimensions of the reaction–gene matrix (rxnGeneMat), if present
%       - Presence and completeness of gene–reaction rules (grRules)
%
%   NOTES:
%       This is a diagnostic/debugging function; it does not modify the model.
%

    % --- Print header for clarity ---
    fprintf('\n--- BASIC DIMENSION CHECKS ---\n');

    % --- Basic counts of model components ---
    fprintf('Reactions: %d\n', numel(model.rxns));
    fprintf('Metabolites: %d\n', numel(model.mets));
    fprintf('Genes: %d\n', numel(model.genes));

    % --- Stoichiometric matrix (S) checks ---
    [mS, nS] = size(model.S);
    if mS ~= numel(model.mets)
        fprintf('Mismatch: S has %d rows, mets = %d\n', mS, numel(model.mets));
    end
    if nS ~= numel(model.rxns)
        fprintf('Mismatch: S has %d cols, rxns = %d\n', nS, numel(model.rxns));
    end

    % --- Reaction–gene matrix checks ---
    if isfield(model, 'rxnGeneMat')
        [mG, nG] = size(model.rxnGeneMat);
        if mG ~= numel(model.rxns)
            fprintf('Mismatch: rxnGeneMat rows = %d, rxns = %d\n', mG, numel(model.rxns));
        end
        if nG ~= numel(model.genes)
            fprintf('Mismatch: rxnGeneMat cols = %d, genes = %d\n', nG, numel(model.genes));
        end
    else
        fprintf('No rxnGeneMat found.\n');
    end

    % --- Gene–reaction rule checks ---
    if isfield(model, 'grRules')
        % Find reactions with missing (empty) gene–reaction rules
        gMissing = find(cellfun(@isempty, model.grRules));
        if ~isempty(gMissing)
            fprintf('Empty grRules in %d reactions\n', numel(gMissing));
        end
    else
        fprintf('No grRules found.\n');
    end
end
