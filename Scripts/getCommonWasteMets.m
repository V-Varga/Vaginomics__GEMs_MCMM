function wasteMets = getCommonWasteMets(model)
% getCommonWasteMets  Return recommended KEGG metabolite IDs for waste exchanges
%
%   wasteMets = getCommonWasteMets(model)
%
%   This function provides a general, organism-independent set of common
%   metabolic byproducts ("waste metabolites") that are often secreted
%   and for which it is sensible to include exchange reactions.
%
%   Inputs:
%     model     COBRA/RAVEN model structure (optional; only used to filter
%               the list to metabolites actually present in the model.mets)
%
%   Outputs:
%     wasteMets Cell array of KEGG IDs (as strings)
%
%   Notes:
%     - The returned list is based on KEGG compound IDs (Cx…).
%     - If model is provided, the list is intersected with model.mets.
%     - Typical entries include CO2, H2O, protons, acetate, lactate, etc.
%
% Example:
%   wasteMets = getCommonWasteMets(LcrispatusBiomassObjDraftModel);
%   allExMets = unique([biomassMets; essentialMets; wasteMets]);

    % Default KEGG IDs for general waste products:
    defaultWaste = { ...
        'C00011', ... % CO2
        'C00001', ... % H2O
        'C00080', ... % H+
        'C00033', ... % Acetate
        'C00186', ... % Lactate
        'C00469', ... % Ethanol
        'C00058', ... % Formate
        'C00022', ... % Pyruvate
        'C00042', ... % Succinate
        'C00154', ... % Propionate
        'C00246', ... % Butyrate
        'C00116', ... % Glycerol
        'C00084', ... % Acetaldehyde
        };

    if nargin < 1 || isempty(model)
        % Return the full default set
        wasteMets = defaultWaste;
    else
        % Intersect with metabolites actually present in the model
        wasteMets = intersect(defaultWaste, model.mets);
    end
end
