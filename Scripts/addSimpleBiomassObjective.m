function model = addSimpleBiomassObjective(model, varargin)
% addSimpleBiomassObjective(model; opts...)
%   Adds a simple biomass objective function to a model.
%
%   Specifically:
%       - Adds a biomass metabolite: biomass_c (in compartment 'c')
%       - Adds a biomass reaction  : BIOMASS_ATP
%         (ATP + H2O -> ADP + Pi + H+ + biomass)
%
%   Options:
%       'ATPperBiomass', <number>   default: 50
%           ATP molecules hydrolyzed per biomass unit
%
%       'Compartment',   <char>     default: 'c'
%           Compartment ID in which biomass metabolite is created
%
%   Input:
%       model    - COBRA/RAVEN-style model structure
%       varargin - optional parameters (see above)
%
%   Output:
%       model    - model with added biomass metabolite and reaction,
%                  and with objective set to the biomass reaction

    %--- Parse optional arguments
    p = inputParser;
    addParameter(p,'ATPperBiomass',50,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(p,'Compartment','c',@(s)ischar(s)||isstring(s));
    parse(p,varargin{:});
    ATPper = p.Results.ATPperBiomass;     % number of ATP per biomass unit
    comp   = char(p.Results.Compartment); % compartment ID

    %--- Ensure the compartment exists in the model
    if ~isfield(model,'comps') || isempty(model.comps)
        model.comps = {comp};
        model.compNames = {comp};
    elseif ~ismember(comp, model.comps)
        model.comps{end+1}    = comp;
        model.compNames{end+1}= comp;
    end
    cIdx = find(strcmp(model.comps, comp), 1);

    %--- Nested helper function: find or create a metabolite
    function metID = ensureMet(cands, niceName)
        % cands: list of possible IDs/names to search for
        % niceName: human-readable metabolite name

        if ~iscell(cands), cands = {cands}; end
        metID = '';

        % Search through candidate IDs
        for k = 1:numel(cands)
            hit = find(strcmpi(model.mets, cands{k}), 1);
            if ~isempty(hit)
                metID = model.mets{hit};
                break;
            end
        end

        % If not found, create a new metabolite entry
        if isempty(metID)
            metID = cands{1};  % use first candidate as ID
            if ~isfield(model,'mets')||isempty(model.mets)
                % Initialize metabolite arrays
                model.mets = {metID};
                model.metNames = {niceName};
                model.metComps = cIdx;
                model.metFormulas = {''};
                model.metCharges  = 0;
            else
                % Append new metabolite
                model.mets{end+1}       = metID;
                model.metNames{end+1}   = niceName;
                model.metComps(end+1,1) = cIdx;

                % Ensure metFormulas and metCharges are sized correctly
                if ~isfield(model,'metFormulas') || numel(model.metFormulas)~=numel(model.mets)-1
                    model.metFormulas = padCell(model.metFormulas, numel(model.mets)-1);
                end
                if ~isfield(model,'metCharges') || numel(model.metCharges)~=numel(model.mets)-1
                    model.metCharges  = padNum(model.metCharges,  numel(model.mets)-1);
                end
                model.metFormulas{end+1}= '';
                model.metCharges(end+1,1)= 0;
            end

            % Expand stoichiometric matrix (S) with new metabolite row
            if ~isfield(model,'S') || isempty(model.S)
                model.S = sparse(numel(model.mets), max(1,numel(model.rxns)));
            else
                model.S(end+1,:) = sparse(1, size(model.S,2));
            end

            % Expand RHS vector (b) with new row
            if ~isfield(model,'b') || isempty(model.b)
                model.b = zeros(numel(model.mets),1);
            else
                model.b(end+1,1) = 0;
            end
        end
    end

    %--- Helper functions to pad cell/num arrays
    function c = padCell(c, n)
        if isempty(c), c = cell(n,1); c(:)={''}; return; end
        c = c(:);
        if numel(c)<n, c(end+1:n,1) = {''}; end
    end
    function v = padNum(v, n)
        if isempty(v), v = zeros(n,1); return; end
        v = v(:);
        if numel(v)<n, v(end+1:n,1) = 0; end
    end

    %--- Standard metabolites: ensure ATP hydrolysis participants exist
    ATP = ensureMet({['atp_' comp], 'ATP', ['ATP[' comp ']'], 'C00002', 'cpd00002'}, 'ATP');
    H2O = ensureMet({['h2o_' comp], 'H2O', ['H2O[' comp ']'], 'C00001', 'water', 'cpd00001'}, 'Water');
    ADP = ensureMet({['adp_' comp], 'ADP', ['ADP[' comp ']'], 'C00008', 'cpd00008'}, 'ADP');
    PI  = ensureMet({['pi_'  comp], 'Pi',  ['PI['  comp ']'], 'C00009', 'phosphate', 'cpd00009'}, 'Phosphate');
    H   = ensureMet({['h_'   comp], 'H+',  ['H['   comp ']'], 'C00080', 'cpd00067'}, 'Proton');

    %--- Biomass metabolite
    biomassID = ['biomass_' comp];
    if ~ismember(biomassID, model.mets)
        ensureMet(biomassID, 'Biomass (aggregate)');
    end

    %--- Biomass reaction setup
    rxnID = 'BIOMASS_ATP';
    if ~isfield(model,'rxns') || isempty(model.rxns)
        % Initialize reaction arrays
        model.rxns = {};
        model.rxnNames = {};
        model.lb = [];
        model.ub = [];
        model.rev = [];
        model.c = [];
        model.equations = {};
    end
    if ~ismember(rxnID, model.rxns)
        % Expand stoichiometric matrix (S) with new reaction column
        if ~isfield(model,'S') || isempty(model.S)
            model.S = sparse(numel(model.mets),1);
        else
            model.S(:,end+1) = sparse(size(model.S,1),1);
        end

        % Add reaction entry
        model.rxns{end+1,1}     = rxnID;
        model.rxnNames{end+1,1} = 'ATP-coupled biomass synthesis';
        model.lb(end+1,1)       = 0;
        model.ub(end+1,1)       = 1000;
        model.rev(end+1,1)      = 0;
        model.c(end+1,1)        = 0;   % objective set later

        % Add stoichiometry: ATP hydrolysis * ATPper -> biomass
        col = numel(model.rxns);
        [~,iATP] = ismember(ATP, model.mets);
        [~,iH2O] = ismember(H2O, model.mets);
        [~,iADP] = ismember(ADP, model.mets);
        [~,iPI ] = ismember(PI,  model.mets);
        [~,iH  ] = ismember(H,   model.mets);
        [~,iBM ] = ismember(biomassID, model.mets);

        model.S(iATP,col) = model.S(iATP,col) - ATPper;
        model.S(iH2O,col) = model.S(iH2O,col) - ATPper;
        model.S(iADP,col) = model.S(iADP,col) + ATPper;
        model.S(iPI ,col) = model.S(iPI ,col) + ATPper;
        model.S(iH  ,col) = model.S(iH  ,col) + ATPper;
        model.S(iBM ,col) = model.S(iBM ,col) + 1;

        % Optional human-readable reaction equation
        model.equations{end+1,1} = sprintf('%g %s + %g %s => %g %s + %g %s + %g %s + %s', ...
            ATPper, ATP, ATPper, H2O, ATPper, ADP, ATPper, PI, ATPper, H, biomassID);
    end

    %--- Set objective to biomass reaction
    obj = zeros(numel(model.rxns),1);
    obj(strcmp(model.rxns, rxnID)) = 1;
    model = setParam(model,'obj',model.rxns,obj);   % RAVEN objective assignment

    %--- Ensure RHS vector (b) is correctly sized
    if ~isfield(model,'b') || numel(model.b)~=numel(model.mets)
        model.b = zeros(numel(model.mets),1);
    end
end
