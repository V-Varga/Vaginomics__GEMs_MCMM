function [model, biomassRxnID] = addPseudoBiomass_trichomonas(model, biomassRxnID, precursors, debugMode, debugFile, csvFile)
% addPseudoBiomass_trichomonas
%
% PURPOSE:
%   Adds a "pseudo biomass" reaction tailored for *Trichomonas vaginalis*.
%   Behaves exactly like the Gram+/Gram- versions, with the same call
%   structure and options.
%
% FEATURES:
%   - Uses partial (case-insensitive) matches on reaction NAMES
%   - Default precursor list reflects known auxotrophies and physiology
%   - Logs matches, reactions, and metabolites
%   - Optional debug logging and CSV export
%   - Stoichiometry = 1 for each metabolite occurrence
%
% USAGE EXAMPLE:
%   [model, bioID] = addPseudoBiomass_trichomonas(model, 'BIOMASS_trich', [], true, ...
%       'Tvaginalis_biomass_debug.txt', 'Tvaginalis_biomass_mets.csv');
%
% INPUTS:
%   model         - GEM in RAVEN/COBRA format
%   biomassRxnID  - String ID for new biomass reaction (default = 'BIOMASS_trich')
%   precursors    - Cell array of precursor keywords (optional; defaults to T. vaginalis list)
%   debugMode     - Boolean; if true, print debug info to console (default = false)
%   debugFile     - Optional path for debug log file ('' = no file)
%   csvFile       - Optional path for CSV export ('' = no file)
%
% OUTPUTS:
%   model         - Updated GEM with biomass reaction added
%   biomassRxnID  - ID of the biomass reaction
%
% NOTES:
%   - Based on literature: Yarlett & Martinez 1986, Müller 1993, 
%     Carlton et al. 2007, Hirt 2013
%
% -------------------------------------------------------------------------

    %% 1. Handle optional inputs
    if nargin < 2 || isempty(biomassRxnID)
        biomassRxnID = 'BIOMASS_trich';
    end
    if nargin < 3 || isempty(precursors)
        % Default precursor list tailored to T. vaginalis
        precursors = {'L-arginine','L-cysteine','L-methionine','L-serine', ...
                      'L-leucine','L-isoleucine','L-valine','L-histidine','L-threonine','L-lysine', ...
                      'L-phenylalanine','L-tyrosine','L-tryptophan','L-glutamine','L-glutamate', ...
                      'L-asparagine','L-aspartate','glycine','ATP','GTP','CTP','UTP','NAD','NADP','FAD','CoA', ...
                      'cholesterol','ergosterol','sterol','phosphatidylcholine','phosphatidylethanolamine', ...
                      'phosphatidylserine','sphingolipid','heme','iron','ferredoxin','thiamine','riboflavin','biotin', ...
                      'glucose','pyruvate','malate'};
    end
    if nargin < 4, debugMode = false; end
    if nargin < 5, debugFile = ''; end
    if nargin < 6, csvFile = ''; end

    %% 2. Call the generic biomass builder
    [model, biomassRxnID] = addPseudoBiomass_debug_all( ...
        model, biomassRxnID, precursors, debugMode, debugFile, csvFile);

end
