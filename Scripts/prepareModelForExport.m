function model = prepareModelForExport(model, outFile)
% prepareModelForExport  Normalize RAVEN model fields for SBML export
%                        Ensures all reaction-, gene-, and rule-level
%                        fields have correct length and type, so that
%                        exportModel() does not fail.
%
%   model = prepareModelForExport(model, outFile)
%
%   Inputs:
%       model   - RAVEN model struct
%       outFile - (optional) filename for SBML export
%
%   Output:
%       model   - updated model with padded/truncated fields
%
%   Description:
%       RAVEN and COBRA models often contain fields of inconsistent length 
%       (e.g. rxnNotes shorter than the number of reactions, or grRules 
%       missing completely). SBML export requires each of these fields to 
%       have exactly one entry per reaction (same length as model.rxns). 
%       This function:
%         • Normalizes reaction-level fields (notes, references, etc.)
%         • Ensures grRules exists and is consistent
%         • Fixes rxnGeneMat so its dimensions match model.rxns and model.genes
%         • Optionally exports the cleaned model to SBML via exportModel()
%
%   See also: exportModel, RAVEN Toolbox

    nRxns = numel(model.rxns); % number of reactions in the model

    %--------------------------------------------------------------
    % Reaction-level fields (rxn-prefixed, eccodes, subSystems)
    %--------------------------------------------------------------
    fnames = fieldnames(model);
    rxnFields = fnames( startsWith(fnames,'rxn') | ...
                        ismember(fnames,{'eccodes','subSystems'}) );

    for i = 1:numel(rxnFields)
        fname = rxnFields{i};
        val   = model.(fname);

        % --------- Handle cell arrays ---------
        if iscell(val)
            if strcmp(fname,'rxnMiriams')
                padVal = {}; % rxnMiriams = cell array of structs
            elseif strcmp(fname,'subSystems')
                % Each subsystem entry must itself be a cell array of strings
                val = cellfun(@(x) forceCellstr(x), val(:), 'UniformOutput', false);
                padVal = {''};
            else
                % General case (rxnReferences, rxnNotes, eccodes, etc.)
                val = cellfun(@char, val(:), 'UniformOutput', false);
                padVal = '';
            end

            % Pad or truncate to match number of reactions
            if numel(val) ~= nRxns
                val = padOrTruncateCell(val,nRxns,padVal);
            end
            model.(fname) = val;

        % --------- Handle numeric/logical arrays ---------
        elseif isnumeric(val) || islogical(val)
            % Ensure numeric arrays are padded with NaN
            model.(fname) = padOrTruncate(val,nRxns,NaN);

        % --------- Unexpected type ---------
        else
            % Fallback: force into cellstr and pad
            warning('Field %s had unexpected type. Forcing into cellstr.',fname);
            tmp = cellstr(string(val(:)));
            model.(fname) = padOrTruncateCell(tmp,nRxns,'');
        end
    end

    %--------------------------------------------------------------
    % grRules field (special case for exportModel)
    %--------------------------------------------------------------
    if ~isfield(model,'grRules')
        % If missing, create empty rules for each reaction
        model.grRules = repmat({''},nRxns,1);
    else
        % Otherwise pad or truncate to correct length
        model.grRules = padOrTruncateCell(model.grRules,nRxns,'');
    end

    %--------------------------------------------------------------
    % Gene-level consistency check
    %--------------------------------------------------------------
    if isfield(model,'genes') && isfield(model,'rxnGeneMat')
        nGenes = numel(model.genes);
        [nRxnMat, nGeneMat] = size(model.rxnGeneMat);

        % Fix mismatch in reaction dimension
        if nRxnMat < nRxns
            model.rxnGeneMat(nRxns,nGeneMat) = 0; % pad with zeros
        elseif nRxnMat > nRxns
            model.rxnGeneMat = model.rxnGeneMat(1:nRxns,:);
        end

        % Fix mismatch in gene dimension
        if nGeneMat < nGenes
            model.rxnGeneMat(:,nGenes) = 0; % pad with zeros
        elseif nGeneMat > nGenes
            model.rxnGeneMat = model.rxnGeneMat(:,1:nGenes);
        end
    end

    %--------------------------------------------------------------
    % Export if requested
    %--------------------------------------------------------------
    if nargin>1 && ~isempty(outFile)
        exportModel(model,outFile);
    end
end

%------------------ Helpers ------------------%
function out = padOrTruncate(vec, n, padVal)
    % Ensure numeric/logical vector has exactly n entries
    if numel(vec) < n
        out = [vec(:); repmat(padVal,n-numel(vec),1)];
    else
        out = vec(1:n);
    end
end

function out = padOrTruncateCell(cellvec, n, padVal)
    % Ensure cell array has exactly n entries
    if numel(cellvec) < n
        out = [cellvec(:); repmat({padVal},n-numel(cellvec),1)];
    else
        out = cellvec(1:n);
    end
end

function out = forceCellstr(x)
    % Convert input to a cell array of strings (always safe for SBML)
    if isempty(x)
        out = {''};
    elseif ischar(x)
        out = {x};
    elseif isstring(x)
        out = cellstr(x);
    elseif iscell(x)
        out = cellfun(@char, x(:), 'UniformOutput', false);
    else
        out = {char(string(x))};
    end
end
