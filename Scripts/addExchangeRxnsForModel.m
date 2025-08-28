function model = addExchangeRxnsForModel(model, metsToAdd)
% addExchangeRxnsForModel  Add exchange reactions to a RAVEN/COBRA model
%
%   model = addExchangeRxnsForModel(model, metsToAdd)
%
%   PURPOSE:
%       Adds exchange reactions for selected metabolites to a GEM.
%       Each reaction allows uptake and secretion with bounds [-1000, 1000].
%
%   INPUTS:
%       model       - struct, GEM in RAVEN or COBRA style
%       metsToAdd   - cell array of metabolite IDs (e.g., {'C00031','C00007'})
%                     Optional. If empty or omitted, exchanges are added for all metabolites.
%
%   OUTPUT:
%       model       - Updated GEM with added exchange reactions
%
%   NOTES:
%       Exchange reactions are named 'EX_<met>' and follow the convention:
%           negative stoichiometry for uptake, positive for secretion.

    % --- Step 1: Default to all metabolites if none specified ---
    if nargin < 2 || isempty(metsToAdd)
        metsToAdd = model.mets; % fallback: all metabolites
        fprintf('No specific list given. Adding exchange for ALL %d metabolites.\n', numel(metsToAdd));
    end

    % --- Step 2: Loop through each metabolite ---
    for i = 1:numel(metsToAdd)
        metID = metsToAdd{i};
        metIdx = find(strcmp(model.mets, metID));
        if isempty(metIdx)
            warning('Metabolite %s not found in model, skipping.', metID);
            continue;
        end

        % --- Step 3: Define exchange reaction ---
        rxnID   = ['EX_' metID];
        rxnName = ['Exchange for ' metID];
        if isfield(model,'rxns') && any(strcmp(model.rxns, rxnID))
            % Skip if exchange already exists
            continue;
        end

        % --- Step 4: Build stoichiometry column ---
        Scol = sparse(numel(model.mets), 1);
        Scol(metIdx) = -1; % convention: negative for uptake

        % --- Step 5: Append reaction to model ---
        model.S(:, end+1)       = Scol;
        model.rxns{end+1,1}     = rxnID;
        model.rxnNames{end+1,1} = rxnName;

        % --- Step 6: Set reaction bounds and objective ---
        model.lb(end+1,1) = -1000;
        model.ub(end+1,1) = 1000;
        model.rev(end+1,1) = 1;
        model.c(end+1,1)  = 0;

        % --- Step 7: Pad optional fields if they exist ---
        if isfield(model,'subSystems');       model.subSystems{end+1,1} = ''; end
        if isfield(model,'eccodes');          model.eccodes{end+1,1} = ''; end
        if isfield(model,'rxnMiriams');       model.rxnMiriams{end+1,1} = struct; end
        if isfield(model,'rxnNotes');         model.rxnNotes{end+1,1} = ''; end
        if isfield(model,'grRules');          model.grRules{end+1,1} = ''; end
        if isfield(model,'rules');            model.rules{end+1,1} = ''; end
        if isfield(model,'equations');        model.equations{end+1,1} = ''; end
        if isfield(model,'rxnReferences');    model.rxnReferences{end+1,1} = ''; end
        if isfield(model,'rxnConfidenceScores'); 
            model.rxnConfidenceScores(end+1,1) = 0; 
        end
        if isfield(model,'rxnFrom');          model.rxnFrom{end+1,1} = ''; end
    end
end
