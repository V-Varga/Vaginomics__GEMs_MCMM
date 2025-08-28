function [model, biomassRxnID] = addPseudoBiomass_debug_all(model, biomassRxnID, precursors, debugMode, debugFile, csvFile)
% addPseudoBiomass_debug_all
%
% PURPOSE:
%   Adds a "pseudo biomass" reaction by collecting ALL metabolites from
%   reactions whose names partially match a given set of precursor names.
%   Useful for exploratory cases when a biomass reaction is missing.
%
% FEATURES:
%   - Uses partial (case-insensitive) matches on reaction NAMES
%   - Includes ALL matches for a precursor (duplicates allowed)
%   - Logs every match (reaction IDs, names, and metabolites)
%   - Optional debug logging to file
%   - Optional CSV export of biomass metabolite list + sources
%   - Stoichiometry = 1 for each occurrence
%
% USAGE EXAMPLE:
%   [model, bioID] = addPseudoBiomass_debug_all(model, 'BIOMASS_ps', [], true, ...
%       'biomass_debug.txt', 'biomass_mets.csv');
%
% INPUTS:
%   model         - GEM in RAVEN/COBRA format
%   biomassRxnID  - String ID for new biomass reaction (default = 'BIOMASS_ps')
%   precursors    - Cell array of precursor keywords (optional; defaults to Gram-positive list)
%   debugMode     - Boolean; if true, print debug info to console (default = false)
%   debugFile     - Optional path for debug log file ('' = no file)
%   csvFile       - Optional path for CSV export ('' = no file)
%
% OUTPUTS:
%   model         - Updated GEM with biomass reaction added
%   biomassRxnID  - ID of the biomass reaction
%
% LIMITATIONS:
%   - Stoichiometries are all set to 1 per occurrence
%   - Biomass composition is NOT biologically accurate; exploratory use only
%
% -------------------------------------------------------------------------

    %% 1. Handle optional inputs
    if nargin < 2 || isempty(biomassRxnID)
        biomassRxnID = 'BIOMASS_ps'; % Default reaction ID if none supplied
    end
    if nargin < 3 || isempty(precursors)
        % Default precursor list (Gram-positive cell components + cofactors)
        precursors = { ...
            'L-alanine', 'L-arginine', 'L-aspartate', 'L-cysteine', ...
            'L-glutamate', 'L-glutamine', 'L-glycine', 'L-histidine', ...
            'L-isoleucine', 'L-leucine', 'L-lysine', 'L-methionine', ...
            'L-phenylalanine', 'L-proline', 'L-serine', 'L-threonine', ...
            'L-tryptophan', 'L-tyrosine', 'L-valine', ...
            'ATP', 'GTP', 'CTP', 'UTP', ...
            'NAD', 'NADP', 'FAD', 'CoA', 'Heme', ...
            'Peptidoglycan', 'Lipoteichoic acid', ...
            'Phosphatidylglycerol', 'Cardiolipin' ...
        };
    end
    if nargin < 4
        debugMode = false; % Default: do not print debug info
    end
    if nargin < 5
        debugFile = ''; % Default: no log file
    end
    if nargin < 6
        csvFile = ''; % Default: no CSV export
    end

    %% 2. Prepare debug logging
    logFID = []; % File identifier (empty unless opened)
    if ~isempty(debugFile)
        logFID = fopen(debugFile, 'w'); % Try opening file for writing
        if logFID == -1
            warning('Could not open debug file "%s". Falling back to console only.', debugFile);
            logFID = []; % Reset to empty if open failed
        end
    end

    % Inline helper function: log to console and/or file
    logMsg = @(msg) ...
        ( ...
            (debugMode && fprintf('%s\n', msg)) || true && ...
            (~isempty(logFID) && fprintf(logFID, '%s\n', msg)) ...
        );

    %% 3. Initialize tracking lists
    metsForBiomass = {};       % Stores metabolite IDs (with duplicates allowed)
    coeffs = [];               % Stores stoichiometric coefficients
    precursorSource = {};      % Records which precursor keyword triggered inclusion
    metNamesList = {};         % Stores metabolite names (for CSV export)

    %% 4. Search for reactions matching each precursor
    for i = 1:numel(precursors)
        target = lower(precursors{i}); % Lowercase for case-insensitive match

        % Find all reactions whose *name* contains this precursor keyword
        matchIdx = find(contains(lower(model.rxnNames), target));

        % Log summary of how many matches found
        if isempty(matchIdx)
            logMsg(sprintf('Precursor "%s" → no matches found', precursors{i}));
            continue;
        else
            logMsg(sprintf('Precursor "%s" → %d matches found', precursors{i}, numel(matchIdx)));
        end

        % Process each matching reaction
        for j = 1:numel(matchIdx)
            rid = model.rxns{matchIdx(j)};      % Reaction ID
            rname = model.rxnNames{matchIdx(j)}; % Reaction name
            logMsg(sprintf('    Match: "%s" (ID: %s)', rname, rid));

            % Extract reactants (negative coefficients in stoichiometric matrix)
            rxnMets = find(model.S(:, matchIdx(j)) < 0);

            % Log and store each metabolite
            for m = 1:numel(rxnMets)
                metID = model.mets{rxnMets(m)};
                metName = model.metNames{rxnMets(m)};
                logMsg(sprintf('        Metabolite: "%s" (ID: %s)', metName, metID));

                % Store information (duplicates kept intentionally)
                metsForBiomass{end+1, 1} = metID;
                metNamesList{end+1, 1} = metName;
                coeffs(end+1, 1) = 1; % Assign stoichiometry = 1
                precursorSource{end+1, 1} = precursors{i};
            end
        end
    end

    %% 5. Add new biomass reaction to model
    % Note: duplicates are preserved; no collapsing of metabolites
    model = addReaction(model, biomassRxnID, ...
        'metaboliteList', metsForBiomass, ...
        'stoichCoeffList', coeffs, ...
        'reversible', false, ...
        'lowerBound', 0, ...
        'upperBound', 1000, ...
        'objectiveCoef', 1);

    %% 6. Export biomass metabolite table to CSV (optional)
    if ~isempty(csvFile)
        T = table(metsForBiomass, metNamesList, coeffs, precursorSource, ...
            'VariableNames', {'Metabolite_ID', 'Metabolite_Name', 'StoichCoeff', 'Source_Precursor'});
        writetable(T, csvFile);
        logMsg(sprintf('CSV table saved to "%s"', csvFile));
    end

    %% 7. Final confirmation message
    logMsg(sprintf('Added biomass reaction "%s" with %d metabolite entries (duplicates kept).', ...
        biomassRxnID, numel(metsForBiomass)));

    %% 8. Close debug file if it was opened
    if ~isempty(logFID)
        fclose(logFID);
    end
end
