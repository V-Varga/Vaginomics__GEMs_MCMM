function model = syncModelStruct(model, logFileName)
% syncModelStruct  
% Synchronize the dimensions of fields in a COBRA/RAVEN model structure.
%
%   model = syncModelStruct(model, logFileName)
%
%   PURPOSE:
%       This function inspects every major field in a COBRA/RAVEN model structure and
%       ensures that its number of rows matches the "true" count for its category:
%           - Reaction-related fields → nRxns (length(model.rxns))
%           - Metabolite-related fields → nMets (length(model.mets))
%           - Gene-related fields → nGenes (length(model.genes))
%       If a field is shorter than expected, it is padded with safe default values.
%       If it is longer, it is trimmed to the correct size.
%       If a field is missing entirely, it is created with safe default values.
%
%       Special handling is provided for:
%           - Numeric vectors (lb, ub, rev, b, metCharges)
%           - MIRIAM annotation fields (rxnMiriams, metMiriams, geneMiriams)
%             → these are cell arrays of structs with 'name' and 'value' fields.
%
%   INPUTS:
%       model        - struct, COBRA/RAVEN model
%       logFileName  - string, optional. Name of the CSV log file where changes are recorded.
%                      If omitted, defaults to 'syncModelStruct_log.csv'.
%
%   OUTPUT:
%       model        - struct, same as input but with all relevant fields padded/truncated/created
%                      so they match their category lengths and have safe defaults.
%
%   LOGGING:
%       All changes are recorded in a CSV file with columns:
%           Field, ExpectedLength, OriginalLength, Action, Details

    % --- Default the log file name if none provided ---
    if nargin < 2
        logFileName = 'syncModelStruct_log.csv';
    end

    % --- Storage for log entries (to be written later) ---
    % Each row will be {fieldName, expectedLen, originalLen, action, details}
    logData = {};

    % --- Determine "truth" lengths for categories ---
    nRxns  = numel(model.rxns);   % number of reactions in the model
    nMets  = numel(model.mets);   % number of metabolites
    nGenes = numel(model.genes);  % number of genes

    % --- Define checklist of fields to inspect ---
    % For each field: name, expected length, and category
    checkList = {
        % Reaction-related fields
        'rxnNames',      nRxns,  'rxn';
        'equations',     nRxns,  'rxn';
        'subSystems',    nRxns,  'rxn';
        'eccodes',       nRxns,  'rxn';
        'rxnReferences', nRxns,  'rxn';
        'rxnConfidenceScores', nRxns, 'rxn';
        'rxnNotes',      nRxns,  'rxn';
        'rxnMiriams',    nRxns,  'rxn'; % MIRIAM: cell array of structs
        'lb',            nRxns,  'rxn'; % numeric lower bounds
        'ub',            nRxns,  'rxn'; % numeric upper bounds
        'rev',           nRxns,  'rxn'; % numeric reversibility flags

        % Metabolite-related fields
        'metNames',      nMets,  'met';
        'b',             nMets,  'met'; % RHS of S·v = b
        'metCharges',    nMets,  'met'; % metabolite charges
        'metFormulas',   nMets,  'met';
        'metMiriams',    nMets,  'met'; % MIRIAM annotations

        % Gene-related fields
        'geneNames',     nGenes, 'gene';
        'geneShortNames',nGenes, 'gene';
        'geneMiriams',   nGenes, 'gene'; % MIRIAM annotations
    };

    % --- Iterate through all fields in checklist ---
    for i = 1:size(checkList,1)
        fName       = checkList{i,1}; % field name (string)
        expectedLen = checkList{i,2}; % required number of entries
        category    = checkList{i,3}; % rxn/met/gene category

        if isfield(model, fName)
            % FIELD ALREADY EXISTS
            currentLen = size(model.(fName), 1);

            % === CASE 1: MIRIAM annotation fields ===
            if ismember(fName, {'rxnMiriams','metMiriams','geneMiriams'})
                emptyStruct = struct('name', {}, 'value', {}); % define "empty" annotation

                if currentLen < expectedLen
                    % Too short → pad with empty structs
                    padSize = expectedLen - currentLen;
                    model.(fName)(end+1:expectedLen, 1) = {emptyStruct};
                    logData(end+1,:) = {fName, expectedLen, currentLen, 'Padded', ...
                        sprintf('+%d empty structs', padSize)};
                elseif currentLen > expectedLen
                    % Too long → trim down
                    model.(fName) = model.(fName)(1:expectedLen, 1);
                    logData(end+1,:) = {fName, expectedLen, currentLen, 'Trimmed', ...
                        sprintf('Removed %d entries', currentLen - expectedLen)};
                end
                continue % Skip further checks (already handled)
            end

            % === CASE 2: Numeric vector fields ===
            if ismember(fName, {'lb','ub','rev','b','metCharges'})
                if currentLen < expectedLen
                    padSize = expectedLen - currentLen;
                    % Assign safe defaults depending on field type
                    switch fName
                        case 'lb', defaultVal = 0;     % lower bound defaults to 0
                        case 'ub', defaultVal = 1000;  % upper bound defaults high
                        case 'rev', defaultVal = 0;    % irreversible by default
                        otherwise, defaultVal = NaN;   % unknown numeric → NaN
                    end
                    model.(fName)(end+1:expectedLen, 1) = defaultVal;
                    logData(end+1,:) = {fName, expectedLen, currentLen, 'Padded', ...
                        sprintf('+%d numeric entries (%g)', padSize, defaultVal)};
                elseif currentLen > expectedLen
                    % Too long → trim
                    model.(fName) = model.(fName)(1:expectedLen, 1);
                    logData(end+1,:) = {fName, expectedLen, currentLen, 'Trimmed', ...
                        sprintf('Removed %d entries', currentLen - expectedLen)};
                end
                continue
            end

            % === CASE 3: Generic cell or numeric fields ===
            if currentLen < expectedLen
                padSize = expectedLen - currentLen;
                if iscell(model.(fName))
                    % For cell arrays: pad with empty strings
                    model.(fName)(end+1:expectedLen, 1) = {''};
                else
                    % For numeric arrays: pad with NaN
                    model.(fName)(end+1:expectedLen, 1) = nan;
                end
                logData(end+1,:) = {fName, expectedLen, currentLen, 'Padded', ...
                    sprintf('+%d entries', padSize)};
            elseif currentLen > expectedLen
                % Too long → trim down
                model.(fName) = model.(fName)(1:expectedLen, 1);
                logData(end+1,:) = {fName, expectedLen, currentLen, 'Trimmed', ...
                    sprintf('Removed %d entries', currentLen - expectedLen)};
            end

        else
            % FIELD IS MISSING COMPLETELY → CREATE with defaults

            % MIRIAM annotations: fill with empty structs
            if ismember(fName, {'rxnMiriams','metMiriams','geneMiriams'})
                emptyStruct = struct('name', {}, 'value', {});
                model.(fName) = repmat({emptyStruct}, expectedLen, 1);
                logData(end+1,:) = {fName, expectedLen, 0, 'Created', 'All empty structs'};

            % Numeric fields with safe defaults
            elseif ismember(fName, {'lb','ub','rev','b','metCharges'})
                switch fName
                    case 'lb', defaultVal = 0;
                    case 'ub', defaultVal = 1000;
                    case 'rev', defaultVal = 0;
                    otherwise, defaultVal = NaN;
                end
                model.(fName) = repmat(defaultVal, expectedLen, 1);
                logData(end+1,:) = {fName, expectedLen, 0, 'Created', ...
                    sprintf('All numeric entries (%g)', defaultVal)};

            % Generic fields: fill with empty strings
            else
                model.(fName) = repmat({''}, expectedLen, 1);
                logData(end+1,:) = {fName, expectedLen, 0, 'Created', 'All empty cells'};
            end
        end
    end

    % --- Write log data to CSV file ---
    headers = {'Field', 'ExpectedLength', 'OriginalLength', 'Action', 'Details'};
    fid = fopen(logFileName, 'w');
    fprintf(fid, '%s,%s,%s,%s,%s\n', headers{:}); % write header row
    for r = 1:size(logData,1)
        fprintf(fid, '%s,%d,%d,%s,%s\n', ...
            logData{r,1}, logData{r,2}, logData{r,3}, logData{r,4}, logData{r,5});
    end
    fclose(fid);

    % --- Print status message to console ---
    fprintf('syncModelStruct: Log saved to "%s"\n', logFileName);
end
