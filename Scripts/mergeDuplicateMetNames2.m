function model = mergeDuplicateMetNames2(model)
% mergeDuplicateMetNames2  Merge metabolites with duplicate names in the same compartment.
%
%   model = mergeDuplicateMetNames2(model)
%
%   This function detects metabolites that share the same (normalized) name 
%   within the same compartment, and merges them into a single metabolite 
%   entry. Stoichiometric coefficients are combined, and redundant entries 
%   are removed. The function runs iteratively until no duplicates remain.
%
%   INPUT:
%       model   COBRA/RAVEN-style model structure with at least fields:
%                   - mets
%                   - metNames
%                   - metComps
%                   - comps
%                   - S (stoichiometric matrix)
%
%   OUTPUT:
%       model   Updated model with duplicates merged.
%
%   Notes:
%     - Name normalization is case-insensitive, trims whitespace, 
%       replaces variant dash characters with '-', collapses multiple 
%       spaces, and strips non-ASCII characters.
%     - Only metabolites in the same compartment can be merged.
%     - Multiple rounds of merging are applied until no duplicates remain.
%
%   Example:
%       model = mergeDuplicateMetNames2(model);

    %--------------------------------------------------------------
    % 1. Ensure metComps is numeric indices (not strings of comp IDs)
    %--------------------------------------------------------------
    if isfield(model,'metComps') && iscell(model.metComps)
        [~, model.metComps] = ismember(model.metComps, model.comps);
    end

    %--------------------------------------------------------------
    % 2. Iteratively normalize names and merge duplicates
    %--------------------------------------------------------------
    changed = true;
    while changed
        %--- Normalize metabolite names
        normNames = lower(strtrim(model.metNames));           % lowercase + trim
        normNames = regexprep(normNames,'[–−]', '-');         % unify dash characters
        normNames = regexprep(normNames,'\s+', ' ');          % collapse multiple spaces
        normNames = regexprep(normNames,'[^ -~]', '');        % strip non-printable ASCII

        %--- Build a unique key: normalized name + compartment
        key = strcat(normNames, "_", string(model.metComps));

        %--- Find duplicate groups
        [~, ~, ic] = unique(key, 'stable');
        changed = false; % reset flag

        %--- Merge duplicates group by group
        for i = 1:max(ic)
            idx = find(ic == i);    % indices of metabolites in this group
            if numel(idx) > 1
                keep = idx(1);      % keep the first occurrence
                dupes = idx(2:end); % list of duplicates to merge/remove

                % Merge stoichiometry (add rows in S for duplicates)
                for d = dupes'
                    model.S(keep,:) = model.S(keep,:) + model.S(d,:);
                end

                % Remove redundant metabolite entries from all relevant fields
                model.S(dupes,:)        = [];
                model.mets(dupes)       = [];
                model.metNames(dupes)   = [];
                if isfield(model,'metFormulas'), model.metFormulas(dupes) = []; end
                if isfield(model,'metCharges'),  model.metCharges(dupes)  = []; end
                if isfield(model,'b'),          model.b(dupes)           = []; end
                if isfield(model,'metComps'),   model.metComps(dupes)    = []; end
                if isfield(model,'metNotes'),   model.metNotes(dupes)    = []; end

                % Since we merged something, another iteration is needed
                changed = true;
            end
        end
    end
end
