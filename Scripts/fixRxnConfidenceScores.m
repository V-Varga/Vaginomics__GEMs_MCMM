function model = fixRxnConfidenceScores(model)
% fixRxnConfidenceScores  Ensure rxnConfidenceScores field is valid
%
%   model = fixRxnConfidenceScores(model)
%
%   PURPOSE:
%       This function checks and corrects the `rxnConfidenceScores` field
%       in a RAVEN/COBRA-style GEM. Confidence scores should be a numeric
%       column vector with one entry per reaction.
%
%       - If the field does not exist → initializes with zeros.
%       - If the field exists but is a cell array → attempts to convert
%         numeric strings or numeric entries into doubles, empty cells → 0.
%       - If the field exists but has wrong size → resets to zeros.
%       - Always outputs a numeric double column vector of size (nRxns × 1).
%
%   INPUT:
%       model   struct   GEM with reactions in `model.rxns`
%
%   OUTPUT:
%       model   struct   GEM with corrected `rxnConfidenceScores`
%
%   NOTES:
%       - This ensures compatibility with downstream tools (e.g. COBRApy).
%       - Confidence scores are often integers (0–4), but stored as doubles.

    if ~isfield(model,'rxnConfidenceScores')
        % --- Case 1: field missing → initialize to zeros ---
        model.rxnConfidenceScores = zeros(numel(model.rxns),1);
        fprintf('Initialized rxnConfidenceScores for %d reactions\n', numel(model.rxns));

    else
        % --- Case 2: field exists, check type ---
        scores = model.rxnConfidenceScores;

        if iscell(scores)
            % --- Subcase: stored as cell array, need conversion ---
            tmp = zeros(numel(model.rxns),1);
            for i = 1:numel(scores)
                if isempty(scores{i})
                    % Empty cell → set to 0
                    tmp(i) = 0;
                elseif isnumeric(scores{i})
                    % Numeric cell → cast to double
                    tmp(i) = double(scores{i});
                elseif ischar(scores{i})
                    % String cell → attempt numeric conversion
                    tmp(i) = str2double(scores{i});
                    if isnan(tmp(i)), tmp(i) = 0; end
                else
                    % Unexpected type → default to 0
                    tmp(i) = 0;
                end
            end
            model.rxnConfidenceScores = tmp;

        else
            % --- Subcase: assume numeric, enforce column double ---
            model.rxnConfidenceScores = double(scores(:));

            % Check for size mismatch
            if numel(model.rxnConfidenceScores) ~= numel(model.rxns)
                model.rxnConfidenceScores = zeros(numel(model.rxns),1);
                fprintf('Reset rxnConfidenceScores to zeros due to size mismatch\n');
            end
        end
    end
end
