function [model, biomassRxnID] = addPseudoBiomass_gramNegative(model, biomassRxnID, precursors, debugMode, debugFile, csvFile)
% addPseudoBiomass_gramNegative  Add a pseudo biomass reaction for Gram-negative GEMs
%
%   [model, biomassRxnID] = addPseudoBiomass_gramNegative(model, biomassRxnID, precursors, debugMode, debugFile, csvFile)
%
% PURPOSE:
%   Constructs a "pseudo biomass" reaction by collecting metabolites from
%   reactions whose names partially match a set of precursor keywords. Useful
%   when no biomass reaction is available in a Gram-negative draft GEM.
%
% INPUTS:
%   model        - GEM structure in RAVEN/COBRA format
%   biomassRxnID - (optional) ID for the new biomass reaction
%                  Default: 'BIOMASS_psGN'
%   precursors   - (optional) cell array of precursor keywords
%                  Default: common amino acids, cofactors, and Gram-negative
%                           cell components
%   debugMode    - (optional) boolean; if true, prints debug info
%                  Default: false
%   debugFile    - (optional) path to debug log file ('' = no file)
%   csvFile      - (optional) path to CSV export ('' = no export)
%
% OUTPUTS:
%   model        - GEM updated with the pseudo biomass reaction
%   biomassRxnID - ID of the biomass reaction added
%
% NOTES:
%   - Stoichiometric coefficients are set to 1 for all collected metabolites.
%   - Duplicates are kept; no merging is done.
%   - CSV export lists metabolites, coefficients, and precursor sources.

    %% 1. Handle optional inputs
    if nargin < 2 || isempty(biomassRxnID)
        biomassRxnID = 'BIOMASS_psGN';
    end
    if nargin < 3 || isempty(precursors)
        % Default precursor list: amino acids, nucleotides, cofactors, lipids
        precursors = {'L-alanine','L-arginine','L-aspartate','L-cysteine', ...
                      'L-glutamate','L-glutamine','L-glycine','L-histidine', ...
                      'L-isoleucine','L-leucine','L-lysine','L-methionine', ...
                      'L-phenylalanine','L-proline','L-serine','L-threonine', ...
                      'L-tryptophan','L-tyrosine','L-valine', ...
                      'ATP','GTP','CTP','UTP', ...
                      'NAD','NADP','FAD','CoA','Heme','Thiamine','Biotin', ...
                      'Peptidoglycan','Lipid A','Lipopolysaccharide','Core oligosaccharide', ...
                      'Phosphatidylethanolamine','Phosphatidylglycerol','Cardiolipin'};
    end
    if nargin < 4, debugMode = false; end
    if nargin < 5, debugFile = ''; end
    if nargin < 6, csvFile = ''; end

    %% 2. Prepare debug logging
    logFID = [];
    if ~isempty(debugFile)
        logFID = fopen(debugFile, 'w');
        if logFID == -1
            warning('Could not open debug file "%s". Logging to console only.', debugFile);
            logFID = [];
        end
    end
    logMsg = @(msg) (debugMode && fprintf('%s\n', msg)) || ...
                     (~isempty(logFID) && fprintf(logFID,'%s\n',msg));

    %% 3. Initialize tracking lists
    metsForBiomass = {};     % metabolite IDs
    metNamesList   = {};     % metabolite names
    precursorSource = {};    % which precursor triggered inclusion
    coeffs = [];             % stoichiometric coefficients

    %% 4. Search reactions for each precursor
    for i = 1:numel(precursors)
        target = lower(precursors{i});
        matchIdx = find(contains(lower(model.rxnNames), target));

        if isempty(matchIdx)
            logMsg(sprintf('Precursor "%s" → no matches found', precursors{i}));
            continue;
        else
            logMsg(sprintf('Precursor "%s" → %d matches found', precursors{i}, numel(matchIdx)));
        end

        % Loop through all matching reactions
        for j = 1:numel(matchIdx)
            rid = model.rxns{matchIdx(j)};
            rname = model.rxnNames{matchIdx(j)};
            logMsg(sprintf('    Match: "%s" (ID: %s)', rname, rid));

            rxnMets = find(model.S(:, matchIdx(j)) < 0); % reactants only
            if isempty(rxnMets), continue; end  % skip if no reactants

            % Add metabolites to biomass list
            for m = 1:numel(rxnMets)
                metID = model.mets{rxnMets(m)};
                metName = model.metNames{rxnMets(m)};

                metsForBiomass  = [metsForBiomass; {metID}];
                metNamesList     = [metNamesList; {metName}];
                coeffs           = [coeffs; 1];
                precursorSource  = [precursorSource; {precursors{i}}];

                logMsg(sprintf('        Metabolite: "%s" (ID: %s)', metName, metID));
            end
        end
    end

    %% 5. Add the biomass reaction
    model = addReaction(model, biomassRxnID, ...
        'metaboliteList', metsForBiomass, ...
        'stoichCoeffList', coeffs, ...
        'reversible', false, ...
        'lowerBound', 0, ...
        'upperBound', 1000, ...
        'objectiveCoef', 1);

    %% 6. Optional CSV export
    if ~isempty(csvFile)
        T = table(metsForBiomass, metNamesList, coeffs, precursorSource, ...
            'VariableNames', {'Metabolite_ID','Metabolite_Name','StoichCoeff','Source_Precursor'});
        writetable(T, csvFile);
        logMsg(sprintf('CSV table saved to "%s"', csvFile));
    end

    %% 7. Confirmation message
    logMsg(sprintf('Added biomass reaction "%s" with %d metabolite entries (duplicates kept).', ...
        biomassRxnID, numel(metsForBiomass)));

    %% 8. Close debug file if open
    if ~isempty(logFID)
        fclose(logFID);
    end
end
