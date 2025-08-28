function model = convertModelSEEDjson(template)
% convertModelSEEDjson
%   Converts a ModelSEED JSON template structure into a simplified
%   model structure compatible with COBRA/RAVEN-style models.
%
%   Input:
%       template - structure decoded from a ModelSEED JSON template
%
%   Output:
%       model    - simplified model structure with fields:
%                  .id, .rxns, .rxnNames, .lb, .ub, .rev, .c
%                  .mets, .metNames, .metComps
%                  .comps, .compNames
%                  .S (initialized as sparse placeholder)

    % Model identifier
    model.id       = template.id;   % copy the ID field from template

    % Reactions
    model.rxns     = {template.reactions.id}';        % list of reaction IDs
    model.rxnNames = {template.reactions.name}';      % reaction names
    model.lb       = arrayfun(@(r) r.lower_bound, template.reactions)';  % lower bounds
    model.ub       = arrayfun(@(r) r.upper_bound, template.reactions)';  % upper bounds
    model.rev      = double(model.lb < 0);            % reversible flag if lb < 0
    model.c        = zeros(numel(model.rxns),1);      % objective coefficients (all 0)

    % Metabolites
    model.mets     = {template.compounds.id}';        % metabolite IDs
    model.metNames = {template.compounds.name}';      % metabolite names

    % Compound-to-compartment mapping
    model.compComps = {template.compcompounds.templatecompound_ref}';   % references to compounds
    model.compRefs  = {template.compcompounds.templatecompartment_ref}';% references to compartments

    % For now assign all metabolites to compartment 1 (dummy placeholder)
    model.metComps = ones(numel(model.mets),1);

    % Compartments
    model.comps     = {template.compartments.id}';    % compartment IDs
    model.compNames = {template.compartments.name}';  % compartment names

    % Stoichiometric matrix (initialized as all zeros, sparse)
    model.S = sparse(numel(model.mets), numel(model.rxns));

    % Ensure that all arrays are column vectors
    model.lb  = model.lb(:);
    model.ub  = model.ub(:);
    model.c   = model.c(:);
    model.rev = model.rev(:);
end
