function wasteMets = getCommunityExchangeMets()
% getCommunityExchangeMets  Recommended metabolites for community exchanges
%
%   wasteMets = getCommunityExchangeMets()
%
%   PURPOSE:
%       Returns a curated list of KEGG metabolite IDs that are commonly
%       used in community metabolic modeling (MCMMs) for exchange reactions.
%       Includes fermentation products, nitrogen sources, gases, SCFAs,
%       and amino acids commonly involved in cross-feeding.
%
%   OUTPUT:
%       wasteMets  - Cell array of KEGG metabolite IDs

    wasteMets = {
        % Carbon sources / fermentation products
        'C00186'; % Lactate
        'C00033'; % Acetate
        'C00469'; % Ethanol
        'C00042'; % Succinate

        % Nitrogen sources
        'C00014'; % Ammonium
        'C00086'; % Urea

        % Gases
        'C00007'; % Oxygen
        'C00011'; % Carbon dioxide
        'C00282'; % Hydrogen

        % Short-chain fatty acids (SCFAs) / fermentation byproducts
        'C00058'; % Formate
        'C00163'; % Propionate
        'C00246'; % Butyrate

        % Amino acids (common cross-feeding metabolites)
        'C00041'; % Alanine
        'C00025'; % Glutamate
        'C00049'; % Aspartate
        'C00065'; % Serine
        'C00097'; % Cysteine
    };
end
